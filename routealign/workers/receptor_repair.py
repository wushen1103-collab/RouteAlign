#!/usr/bin/env python3
"""PB-triggered bounded repair for receptor-rescored Core-Graft predictions."""
from __future__ import annotations
from routealign import settings

import argparse
import atexit
import hashlib
import json
import os
import pathlib
import re
import resource
import shutil
import signal
import subprocess
import tempfile
import time

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign

import routealign.workers.receptor_score as rescore
from routealign.workers.repair import pbcheck


ROOT = settings.workspace()
INPUT = ROOT / "runs/receptor_score"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_RECEPTOR_RESCORE_REPAIR_OUTPUT",
    str(ROOT / "runs/receptor_repair"),
)).resolve()
SMINA = ROOT / "tools/smina"
EXPECTED_SMINA_SHA256 = "ffe5e1e78c947f76b0df8805e2c54383d0bbaf2e827a633b643a708cf682a958"
CPU = settings.cpu_index()
TRUST_REGION_A = 1.25
ACTIVE: pathlib.Path | None = None


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: pathlib.Path, payload: dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def clean() -> None:
    global ACTIVE
    if ACTIVE is not None and ACTIVE.exists():
        shutil.rmtree(ACTIVE)
    ACTIVE = None


def stop(sig, _frame) -> None:
    clean()
    raise SystemExit(128 + sig)


atexit.register(clean)
signal.signal(signal.SIGTERM, stop)
signal.signal(signal.SIGINT, stop)


def load_single(path: pathlib.Path) -> Chem.Mol:
    molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if mol is not None]
    if len(molecules) != 1 or molecules[0].GetNumConformers() != 1:
        raise RuntimeError("invalid SDF: " + str(path))
    return molecules[0]


def graph_identity(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def expected(index: int) -> tuple[dict, pathlib.Path, pathlib.Path]:
    generator_spec, _row, seed, protein = rescore.expected(index)
    pair = generator_spec["pair"]
    pose = INPUT / pair / f"seed_{seed}" / "selected.sdf"
    metadata = INPUT / pair / f"seed_{seed}" / "metadata.json"
    for path in (pose, metadata, protein, SMINA):
        path.resolve().relative_to(ROOT)
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(path)
    if sha256(SMINA) != EXPECTED_SMINA_SHA256:
        raise RuntimeError("Smina identity drift")
    spec = {"task_index": index, "pair": pair, "seed": seed,
            "input_pose_sha256": sha256(pose), "generation_metadata_sha256": sha256(metadata),
            "protein_sha256": sha256(protein), "smina_sha256": EXPECTED_SMINA_SHA256,
            "generator_reads_native": False, "repair_reads_native": False,
            "trust_region_A": TRUST_REGION_A, "threads_max": 1, "cpu_affinity": [CPU],
            "cpu_reserve_min": 30, "gpu_hidden": True, "rdkit_version": rdBase.rdkitVersion}
    return spec, pose, protein


def valid(path: pathlib.Path, spec: dict) -> bool:
    if not path.is_file():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return payload.get("status") == "success" and all(payload.get(k) == v for k, v in spec.items())


def run(index: int, audit_only: bool = False, fault_after_workload: bool = False,
        hold_after_workload: float = 0.0) -> dict:
    global ACTIVE
    if pathlib.Path.cwd().resolve() != ROOT or len(os.sched_getaffinity(0)) < 1:
        raise RuntimeError("scope or CPU reserve precondition failed")
    spec, pose, protein = expected(index)
    if audit_only:
        return {**spec, "status": "audit_pass"}
    os.sched_setaffinity(0, [CPU])
    resource.setrlimit(resource.RLIMIT_AS, (6 * 1024**3, 6 * 1024**3))
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or os.getpid() != os.getpgrp():
        raise RuntimeError("GPU or process-group isolation missing")
    task_root = OUTPUT / spec["pair"] / f"seed_{spec['seed']}"
    result_path = task_root / "metadata.json"
    if valid(result_path, spec) and (task_root / "selected.sdf").is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if payload["selected_pose_sha256"] != sha256(task_root / "selected.sdf"):
            raise RuntimeError("resumed selected pose drift")
        return {**payload, "status": "resumed"}
    pre_pb, pre_failed = pbcheck(pose, protein)
    task_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.", dir=task_root.parent))
    ACTIVE = stage
    try:
        selected = stage / "selected.sdf"
        source = load_single(pose)
        movement = 0.0
        affinity = None
        repaired_failed = []
        if all(pre_pb.values()):
            shutil.copy2(pose, selected)
            decision = "preserve_valid"
        else:
            candidate = stage / "candidate.sdf"
            command = [str(SMINA), "-r", str(protein), "-l", str(pose), "--minimize",
                       "--minimize_iters", "500", "--seed", str(spec["seed"]),
                       "--cpu", "1", "-o", str(candidate)]
            completed = subprocess.run(
                command, cwd=stage, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                timeout=240, env={**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1",
                                  "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                                  "NUMEXPR_NUM_THREADS": "1"},
            )
            if completed.returncode != 0 or not candidate.is_file():
                raise RuntimeError("Smina repair failed: " + completed.stdout[-2000:])
            hits = re.findall(r"Affinity:\s*([-+0-9.eE]+)", completed.stdout)
            affinity = None if not hits else float(hits[-1])
            repaired = load_single(candidate)
            if graph_identity(repaired) != graph_identity(source):
                raise RuntimeError("repair changed graph/stereo")
            movement = float(rdMolAlign.CalcRMS(source, repaired, maxMatches=1000000))
            repaired_pb, repaired_failed = pbcheck(candidate, protein)
            if all(repaired_pb.values()) and movement <= TRUST_REGION_A:
                shutil.copy2(candidate, selected)
                decision = "repair_accepted"
            else:
                shutil.copy2(pose, selected)
                decision = "repair_rejected"
        if fault_after_workload:
            raise SystemExit(97)
        if hold_after_workload > 0:
            time.sleep(hold_after_workload)
        final_pb, final_failed = pbcheck(selected, protein)
        payload = {**spec, "status": "success", "decision": decision,
                   "repair_triggered": bool(pre_failed), "pre_failed": pre_failed,
                   "movement_A": movement, "affinity_kcal_mol": affinity,
                   "repaired_failed": repaired_failed, "final_PB_pass": all(final_pb.values()),
                   "final_failed": final_failed, "selected_pose_sha256": sha256(selected), "SOTA": False}
        atomic_json(stage / "metadata.json", payload)
        os.replace(stage, task_root)
        ACTIVE = None
        return payload
    finally:
        if ACTIVE is not None:
            clean()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--fault-after-workload", action="store_true")
    parser.add_argument("--hold-after-workload", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(run(args.index, args.audit_only, args.fault_after_workload,
                         args.hold_after_workload), sort_keys=True))


if __name__ == "__main__":
    main()
