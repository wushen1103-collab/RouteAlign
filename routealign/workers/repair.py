#!/usr/bin/env python3
"""Self-contained PB-triggered bounded Smina repair for independent Core-Graft v5."""
from __future__ import annotations
from routealign import settings

import argparse
import atexit
import csv
import hashlib
import json
import os
import pathlib
import re
import shutil
import signal
import subprocess
import tempfile
import time

from rdkit import Chem, rdBase


ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
INPUT = ROOT / "runs/graft"
DATASET = ROOT / "data/AlignDockBench"
SMINA = ROOT / "tools/smina"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ADB_CROSSDOCK_REPAIR_OUTPUT",
    str(ROOT / "runs/repair"),
)).resolve()
SEEDS = (101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163,
         167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229)
EXPECTED_MANIFEST_SHA256 = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
EXPECTED_SMINA_SHA256 = "ffe5e1e78c947f76b0df8805e2c54383d0bbaf2e827a633b643a708cf682a958"
ACTIVE = None


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def clean():
    global ACTIVE
    if ACTIVE and ACTIVE.exists():
        shutil.rmtree(ACTIVE)
    ACTIVE = None


def stop(sig, _):
    clean()
    raise SystemExit(128 + sig)


atexit.register(clean)
signal.signal(signal.SIGTERM, stop)
signal.signal(signal.SIGINT, stop)


def manifest_rows():
    if sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256 or sha256(SMINA) != EXPECTED_SMINA_SHA256:
        raise RuntimeError("frozen input identity drift")
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["final_primary_eligible"] == "True"]
    rows.sort(key=lambda row: (row["family"], row["query_id"]))
    if len(rows) != 44:
        raise RuntimeError("cohort mismatch")
    return rows


def tasks():
    return [(i, row, seed) for i, (row, seed) in enumerate((r, s) for r in manifest_rows() for s in SEEDS)]


def load_single(path):
    molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if mol is not None]
    if len(molecules) != 1 or molecules[0].GetNumConformers() != 1:
        raise RuntimeError("molecule parse failed: " + str(path))
    return molecules[0]


def graph_identity(mol):
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def pbcheck(pose, protein):
    from posebusters import PoseBusters
    frame = PoseBusters(config="dock", max_workers=0).bust(
        str(pose), mol_cond=str(protein), full_report=False
    )
    if len(frame) != 1:
        raise RuntimeError("unexpected PoseBusters row count")
    checks = {str(key): bool(value) for key, value in frame.iloc[0].to_dict().items()}
    return checks, [key for key, value in checks.items() if not value]


def run(index, fault_after_workload=False, hold_after_workload=0.0, audit_only=False):
    global ACTIVE
    if pathlib.Path.cwd().resolve() != ROOT.resolve() or not 0 <= index < 1100 or (os.cpu_count() or 0) < 1:
        raise RuntimeError("scope/resources mismatch")
    _, row, seed = tasks()[index]
    pair = f"{row['family']}__{row['query_id']}"
    pose = INPUT / pair / f"seed_{seed}" / "pose_pocket_hybrid.sdf"
    metadata = INPUT / pair / f"seed_{seed}" / "metadata.json"
    protein = DATASET / row["family"] / row["template_id"] / "protein.pdb"
    for path in (pose, metadata, protein):
        path.resolve().relative_to(ROOT.resolve())
        if not path.is_file():
            raise FileNotFoundError(path)
    expected = {
        "task_index": index,
        "family": row["family"],
        "query_id": row["query_id"],
        "seed": seed,
        "manifest_sha256": EXPECTED_MANIFEST_SHA256,
        "smina_sha256": EXPECTED_SMINA_SHA256,
        "input_pose_sha256": sha256(pose),
        "generation_metadata_sha256": sha256(metadata),
        "protein_sha256": sha256(protein),
        "worker_sha256": sha256(__file__),
    }
    if audit_only:
        return {**expected, "status": "audit_pass", "threads_max": 1, "gpu_hidden": True}
    task_root = OUTPUT / pair / f"seed_{seed}"
    result_path = task_root / "result.json"
    if result_path.is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if all(payload.get(key) == value for key, value in expected.items()):
            return {**expected, "status": "resumed", "decision": payload["decision"]}
        raise RuntimeError("resume identity drift")

    pre_pb, pre_failed = pbcheck(pose, protein)
    task_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.", dir=task_root.parent))
    ACTIVE = stage
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
                   "--minimize_iters", "500", "--seed", str(seed), "--cpu", "1", "-o", str(candidate)]
        completed = subprocess.run(
            command, cwd=stage, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=240, env={**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1",
                              "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
        )
        if completed.returncode != 0 or not candidate.is_file():
            raise RuntimeError("smina failed: " + completed.stdout[-2000:])
        hits = re.findall(r"Affinity:\s*([-+0-9.eE]+)", completed.stdout)
        affinity = None if not hits else float(hits[-1])
        repaired = load_single(candidate)
        if graph_identity(repaired) != graph_identity(source):
            raise RuntimeError("graph/stereo changed")
        from rdkit.Chem import rdMolAlign
        movement = float(rdMolAlign.CalcRMS(source, repaired, maxMatches=1000000))
        repaired_pb, repaired_failed = pbcheck(candidate, protein)
        if all(repaired_pb.values()) and movement <= 1.0:
            shutil.copy2(candidate, selected)
            decision = "repair_accepted"
        else:
            shutil.copy2(pose, selected)
            decision = "repair_rejected"
    if fault_after_workload:
        clean()
        raise SystemExit(97)
    if hold_after_workload > 0:
        time.sleep(hold_after_workload)
    final_pb, final_failed = pbcheck(selected, protein)
    payload = {
        **expected,
        "status": "success",
        "decision": decision,
        "repair_triggered": bool(pre_failed),
        "pre_failed": pre_failed,
        "movement_A": movement,
        "affinity_kcal_mol": affinity,
        "trust_region_A": 1.0,
        "repaired_failed": repaired_failed,
        "final_PB_pass": all(final_pb.values()),
        "final_failed": final_failed,
        "selected_sha256": sha256(selected),
        "query_native_ligand_read": False,
        "target_pose_metric_read": False,
        "rdkit_version": rdBase.rdkitVersion,
        "SOTA": False,
    }
    atomic_json(stage / "result.json", payload)
    os.replace(stage, task_root)
    ACTIVE = None
    return {**expected, "status": "computed", "decision": decision,
            "movement_A": movement, "final_PB_pass": payload["final_PB_pass"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--fault-after-workload", action="store_true")
    parser.add_argument("--hold-after-workload", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(run(args.index, args.fault_after_workload, args.hold_after_workload, args.audit_only), sort_keys=True))


if __name__ == "__main__":
    main()
