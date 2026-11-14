#!/usr/bin/env python3
"""Receptor-conditioned, outcome-blind rescoring of Core-Graft candidates."""
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
import resource
import shutil
import signal
import subprocess
import tempfile
import time

import numpy as np
from rdkit import Chem, rdBase

import routealign.workers.graft as generator
import routealign.workers.graft_fallback as fast_generator


ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
PREREG = ROOT / "config/protocol.json"
DATASET = ROOT / "data/AlignDockBench"
SMINA = ROOT / "tools/smina"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_RECEPTOR_RESCORE_OUTPUT",
    str(ROOT / "runs/receptor_score"),
)).resolve()
EXPECTED_MANIFEST_SHA256 = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
EXPECTED_SMINA_SHA256 = "ffe5e1e78c947f76b0df8805e2c54383d0bbaf2e827a633b643a708cf682a958"
SEEDS = (101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163,
         167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229)
CPU = settings.cpu_index()
STRAIN_WEIGHT = 0.10
SHAPE_WEIGHT = 0.05
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


def load_specs() -> list[dict]:
    if sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256 or sha256(SMINA) != EXPECTED_SMINA_SHA256:
        raise RuntimeError("frozen dependency identity drift")
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["final_primary_eligible"] == "True"]
    by_pair = {f"{row['family']}__{row['query_id']}": row for row in rows}
    specs = [{**item, "row": by_pair[item["pair"]]}
             for item in prereg["pairs"] if not item["hard_low_2d"]]
    if len(specs) != 18:
        raise RuntimeError("non-low-2D cohort drift")
    return specs


def tasks() -> list[tuple[int, dict, int]]:
    return [(index, item, seed) for index, (item, seed) in enumerate(
        (item, seed) for item in load_specs() for seed in SEEDS
    )]


def coordinate_key(mol: Chem.Mol) -> str:
    xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=np.float64)
    rounded = np.round(xyz, 4).astype("<f8", copy=False)
    return hashlib.sha256(rounded.tobytes()).hexdigest()


def smina_score(pose: pathlib.Path, protein: pathlib.Path, seed: int) -> tuple[float, str]:
    command = [str(SMINA), "-r", str(protein), "-l", str(pose), "--score_only",
               "--cpu", "1", "--seed", str(seed)]
    completed = subprocess.run(
        command, cwd=pose.parent, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=120, env={**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1",
                          "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                          "NUMEXPR_NUM_THREADS": "1"},
    )
    hits = re.findall(r"Affinity:\s*([-+0-9.eE]+)", completed.stdout)
    if completed.returncode != 0 or len(hits) != 1:
        raise RuntimeError({"returncode": completed.returncode, "output": completed.stdout[-2000:]})
    return float(hits[0]), " ".join(command)


def expected(index: int) -> tuple[dict, dict, int, pathlib.Path]:
    table = tasks()
    if not 0 <= index < len(table):
        raise IndexError(index)
    _, item, seed = table[index]
    row = item["row"]
    protein = DATASET / row["family"] / row["template_id"] / "protein.pdb"
    for path in (MANIFEST, PREREG, SMINA, protein):
        path.resolve().relative_to(ROOT)
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(path)
    spec = {
        "task_index": index, "pair": item["pair"], "pair_index": item["pair_index"], "seed": seed,
        "hard_low_2d": False, "hard_flex": bool(item["hard_flex"]),
        "hard_scaffold": bool(item["hard_scaffold"]),
        "manifest_sha256": EXPECTED_MANIFEST_SHA256, "preregister_sha256": sha256(PREREG),
        "generator_sha256": sha256(pathlib.Path(generator.__file__)),
        "smina_sha256": EXPECTED_SMINA_SHA256, "protein_sha256": sha256(protein),
        "selection_rule": "minimum Smina affinity + 0.10*relative intramolecular energy per heavy atom + 0.05*shape distance",
        "strain_weight": STRAIN_WEIGHT, "shape_weight": SHAPE_WEIGHT,
        "generator_reads_native": False, "selector_reads_native": False,
        "template_receptor_read": True, "query_native_ligand_read": False,
        "threads_max": 1, "cpu_affinity": [CPU], "cpu_reserve_min": 30,
        "gpu_hidden": True, "external_call_timeout_seconds": 120,
        "rdkit_version": rdBase.rdkitVersion,
    }
    if item["pair"] == "2ica_2IC_A_1__3m6f":
        spec.update({
            "candidate_schedule": "deterministic_nonrandom_timeout_fallback",
            "fast_generator_sha256": sha256(pathlib.Path(fast_generator.__file__)),
            "coordinate_equivalence_to_full_schedule": False,
            "fallback_reason": "full candidate schedule has a preregistered reproducible runtime timeout",
        })
    return spec, row, seed, protein


def valid(path: pathlib.Path, spec: dict) -> bool:
    if not path.is_file():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return payload.get("status") == "success" and all(payload.get(k) == v for k, v in spec.items())


def run(index: int, audit_only: bool = False, fault_after_scoring: bool = False,
        hold_after_scoring: float = 0.0) -> dict:
    global ACTIVE
    if pathlib.Path.cwd().resolve() != ROOT or len(os.sched_getaffinity(0)) < 1:
        raise RuntimeError("scope or CPU reserve precondition failed")
    spec, row, seed, protein = expected(index)
    if audit_only:
        return {**spec, "status": "audit_pass"}
    os.sched_setaffinity(0, [CPU])
    resource.setrlimit(resource.RLIMIT_AS, (6 * 1024**3, 6 * 1024**3))
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or os.getpid() != os.getpgrp():
        raise RuntimeError("GPU or process-group isolation missing")
    task_root = OUTPUT / spec["pair"] / f"seed_{seed}"
    result_path = task_root / "metadata.json"
    if valid(result_path, spec) and (task_root / "selected.sdf").is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if payload["selected_pose_sha256"] != sha256(task_root / "selected.sdf"):
            raise RuntimeError("resumed pose checksum drift")
        return {**payload, "status": "resumed"}
    task_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.", dir=task_root.parent))
    ACTIVE = stage
    try:
        active_generator = fast_generator if spec["pair"] == "2ica_2IC_A_1__3m6f" else generator
        candidates, _template = active_generator.build_candidates(row, seed)
        active_generator.rank_candidates(candidates)
        minimum_energy = min(float(c["energy_per_heavy_atom"]) for c in candidates)
        unique: dict[str, dict] = {}
        for candidate_index, candidate in enumerate(candidates):
            key = coordinate_key(candidate["mol"])
            record = unique.setdefault(key, {"representative": candidate, "candidate_indices": [],
                                             "mapping_indices": [], "conformer_indices": []})
            record["candidate_indices"].append(candidate_index)
            record["mapping_indices"].append(int(candidate["mapping_index"]))
            record["conformer_indices"].append(int(candidate["conformer_index"]))
        scored = []
        for unique_index, (key, record) in enumerate(sorted(unique.items())):
            candidate = record["representative"]
            pose = stage / f"candidate_{unique_index:03d}.sdf"
            Chem.MolToMolFile(candidate["mol"], str(pose))
            affinity, command = smina_score(pose, protein, seed)
            relative_strain = max(0.0, float(candidate["energy_per_heavy_atom"]) - minimum_energy)
            selection_score = affinity + STRAIN_WEIGHT * relative_strain + SHAPE_WEIGHT * float(candidate["shape_distance"])
            scored.append({"unique_index": unique_index, "coordinate_key": key,
                           "pose_path": pose, "pose_sha256": sha256(pose),
                           "affinity_kcal_mol": affinity, "relative_strain_per_heavy_atom": relative_strain,
                           "shape_distance": float(candidate["shape_distance"]),
                           "selection_score": float(selection_score),
                           "mapping_indices": sorted(set(record["mapping_indices"])),
                           "conformer_indices": sorted(set(record["conformer_indices"])),
                           "candidate_indices": record["candidate_indices"],
                           "smina_command": command})
        selected = min(scored, key=lambda x: (x["selection_score"], x["coordinate_key"]))
        selected_path = stage / "selected.sdf"
        shutil.copy2(selected["pose_path"], selected_path)
        ready = {"task_index": index, "pair": spec["pair"], "seed": seed,
                 "candidate_count": len(candidates), "unique_pose_count": len(scored),
                 "selected_unique_index": selected["unique_index"],
                 "selected_pose_sha256": sha256(selected_path),
                 "cpu_affinity": sorted(os.sched_getaffinity(0)), "gpu_hidden": True}
        atomic_json(OUTPUT / f"ready_{index}_{os.getpid()}.json", ready)
        if fault_after_scoring:
            raise SystemExit(97)
        if hold_after_scoring > 0:
            time.sleep(hold_after_scoring)
        for item in scored:
            item.pop("pose_path")
        payload = {**spec, "status": "success", "candidate_count": len(candidates),
                   "unique_pose_count": len(scored), "minimum_energy_per_heavy_atom": minimum_energy,
                   "scores": scored, "selected_unique_index": selected["unique_index"],
                   "selected_pose_sha256": sha256(selected_path), "SOTA": False}
        atomic_json(stage / "metadata.json", payload)
        if task_root.exists():
            if valid(result_path, spec):
                return {**json.loads(result_path.read_text()), "status": "resumed"}
            raise RuntimeError("existing result identity drift")
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
    parser.add_argument("--fault-after-scoring", action="store_true")
    parser.add_argument("--hold-after-scoring", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(run(args.index, args.audit_only, args.fault_after_scoring,
                         args.hold_after_scoring), sort_keys=True))


if __name__ == "__main__":
    main()
