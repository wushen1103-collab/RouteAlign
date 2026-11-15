#!/usr/bin/env python3
"""Outcome-blind group medoid selector for the frozen full 44-pair expansion."""
from __future__ import annotations
from routealign import settings
from routealign.atom_correspondence import matrix
import routealign.atom_correspondence as correspondence

import argparse
import atexit
import hashlib
import json
import os
import pathlib
import shutil
import signal
import tempfile
import time

import numpy as np
from rdkit import Chem


ROOT = settings.workspace()
INPUT = ROOT / "runs/bounded_repair"
PREREG = ROOT / "config/protocol.json"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ADB_MEDOID_OUTPUT",
    str(ROOT / "runs/medoid_select"),
)).resolve()
ACTIVE = None


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, payload):
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


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


def spec():
    data = json.loads(PREREG.read_text(encoding="utf-8"))
    if not data.get("stagewise_frozen_before_full_expansion_batch") or not data.get("method_fixed_after_low2d_development"):
        raise RuntimeError("preregistration mismatch")
    pairs = data["pairs"]
    groups = [tuple(group) for group in data["outer_seed_groups"]]
    if len(pairs) != 44 or len(groups) != 5 or any(len(group) != 5 for group in groups):
        raise RuntimeError("frozen design mismatch")
    return pairs, groups


def tasks():
    pairs, groups = spec()
    return [(index, pair_spec, group_index, group)
            for index, (pair_spec, group_index, group) in enumerate(
                (p, gi, g) for p in pairs for gi, g in enumerate(groups, start=1)
            )]


def coordinates(path):
    molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if mol is not None]
    if len(molecules) != 1 or molecules[0].GetNumConformers() != 1:
        raise RuntimeError("pose parse failed: " + str(path))
    return np.asarray(molecules[0].GetConformer().GetPositions(), dtype=float)


def run(index, audit_only=False, fault_after_selection=False, hold_after_selection=0.0):
    global ACTIVE
    table = tasks()
    if pathlib.Path.cwd().resolve() != ROOT.resolve() or not 0 <= index < len(table):
        raise RuntimeError("scope/index mismatch")
    _, pair_spec, group_index, seeds = table[index]
    pair = pair_spec["pair"]
    poses = [INPUT / pair / f"seed_{seed}" / "selected.sdf" for seed in seeds]
    repairs = [INPUT / pair / f"seed_{seed}" / "result.json" for seed in seeds]
    for path in (*poses, *repairs):
        path.resolve().relative_to(ROOT.resolve())
        if not path.is_file():
            raise FileNotFoundError(path)
    expected = {
        "task_index": index,
        "pair": pair,
        "group_index": group_index,
        "candidate_seeds": list(seeds),
        "hard_low_2d": bool(pair_spec["hard_low_2d"]),
        "hard_flex": bool(pair_spec["hard_flex"]),
        "hard_scaffold": bool(pair_spec["hard_scaffold"]),
        "preregister_sha256": sha256(PREREG),
        "candidate_pose_sha256": [sha256(path) for path in poses],
        "repair_result_sha256": [sha256(path) for path in repairs],
        "mapping_sha256": sha256(pathlib.Path(correspondence.__file__)),
                "worker_sha256": sha256(__file__),
    }
    if audit_only:
        return {**expected, "status": "audit_pass", "selector_reads_native": False,
                "endpoint_metrics_read": False, "threads_max": 1, "gpu_hidden": True}
    task_root = OUTPUT / pair / f"group_{group_index}"
    result_path = task_root / "result.json"
    if result_path.is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if all(payload.get(key) == value for key, value in expected.items()):
            return {**expected, "status": "resumed", "selected_seed": payload["selected_seed"]}
        raise RuntimeError("resume identity drift")
    molecules = [next(m for m in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if m is not None) for path in poses]
    distance_matrix, mapping_audit = matrix(molecules)
    distances = {(a,b): float(distance_matrix[i,j]) for i,a in enumerate(seeds) for j,b in enumerate(seeds)}
    ranking = []
    for seed in seeds:
        other = [distances[seed, candidate] for candidate in seeds if candidate != seed]
        ranking.append((float(np.median(other)), float(np.mean(other)), seed))
    selected_seed = min(ranking)[2]
    if fault_after_selection:
        raise SystemExit(97)
    if hold_after_selection > 0:
        time.sleep(hold_after_selection)
    task_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.", dir=task_root.parent))
    ACTIVE = stage
    source = INPUT / pair / f"seed_{selected_seed}" / "selected.sdf"
    shutil.copy2(source, stage / "selected.sdf")
    payload = {
        **expected,
        "status": "success",
        "selected_seed": selected_seed,
        "selected_pose_sha256": sha256(stage / "selected.sdf"),
        "chemical_mapping_audit": mapping_audit,
        "distance_definition": "minimum connectivity-preserving full-graph, resonance-aware RMSD; no rigid fit",
        "rank_median_mean_seed": [list(item) for item in sorted(ranking)],
        "selector": "receptor-frame coordinate medoid within a frozen disjoint five-seed group",
        "selector_reads_native": False,
        "query_native_ligand_read": False,
        "target_pose_metric_read": False,
        "SOTA": False,
    }
    atomic_json(stage / "result.json", payload)
    os.replace(stage, task_root)
    ACTIVE = None
    return {**expected, "status": "computed", "selected_seed": selected_seed}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--fault-after-selection", action="store_true")
    parser.add_argument("--hold-after-selection", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(run(args.index, args.audit_only, args.fault_after_selection,
                         args.hold_after_selection), sort_keys=True))


if __name__ == "__main__":
    main()
