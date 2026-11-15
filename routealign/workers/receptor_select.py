#!/usr/bin/env python3
"""PB-valid receptor-score/consensus selector locked on seed group 1."""
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
GENERATION = ROOT / "runs/receptor_score"
REPAIR = ROOT / "runs/receptor_repair"
PREREG = ROOT / "config/protocol.json"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_RESCORE_SELECTOR_V2_OUTPUT",
    str(ROOT / "runs/receptor_select"),
)).resolve()
CONSENSUS_WEIGHT = 5.0
ACTIVE: pathlib.Path | None = None


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: pathlib.Path, payload: dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n"); handle.flush(); os.fsync(handle.fileno())


def clean() -> None:
    global ACTIVE
    if ACTIVE is not None and ACTIVE.exists(): shutil.rmtree(ACTIVE)
    ACTIVE = None


def stop(sig, _frame) -> None:
    clean(); raise SystemExit(128 + sig)


atexit.register(clean)
signal.signal(signal.SIGTERM, stop)
signal.signal(signal.SIGINT, stop)


def tasks() -> list[tuple[int, dict, int, tuple[int, ...]]]:
    data = json.loads(PREREG.read_text(encoding="utf-8"))
    if not data.get("stagewise_frozen_before_full_expansion_batch"):
        raise RuntimeError("preregistration mismatch")
    pairs = [item for item in data["pairs"] if not item["hard_low_2d"]]
    groups = [tuple(group) for group in data["outer_seed_groups"]]
    if len(pairs) != 18 or len(groups) != 5 or any(len(group) != 5 for group in groups):
        raise RuntimeError("frozen non-low2d design mismatch")
    return [(index, pair, gi, group) for index, (pair, gi, group) in enumerate(
        (pair, gi, group) for pair in pairs for gi, group in enumerate(groups, start=1))]


def coordinates(path: pathlib.Path) -> np.ndarray:
    mols = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True)
            if mol is not None]
    if len(mols) != 1 or mols[0].GetNumConformers() != 1:
        raise RuntimeError("pose parse failed: " + str(path))
    return np.asarray(mols[0].GetConformer().GetPositions(), dtype=float)


def zscores(values: list[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    std = float(array.std(ddof=0))
    return np.zeros_like(array) if std <= 1e-12 else (array - array.mean()) / std


def run(index: int, audit_only: bool = False, fault_after_selection: bool = False,
        hold_after_selection: float = 0.0) -> dict:
    global ACTIVE
    table = tasks()
    if pathlib.Path.cwd().resolve() != ROOT or not 0 <= index < len(table):
        raise RuntimeError("scope/index mismatch")
    _, pair_spec, group_index, seeds = table[index]
    pair = pair_spec["pair"]
    poses = [REPAIR / pair / f"seed_{seed}" / "selected.sdf" for seed in seeds]
    repair_meta = [REPAIR / pair / f"seed_{seed}" / "metadata.json" for seed in seeds]
    generation_meta = [GENERATION / pair / f"seed_{seed}" / "metadata.json" for seed in seeds]
    for path in (*poses, *repair_meta, *generation_meta):
        path.resolve().relative_to(ROOT)
        if not path.is_file() or path.is_symlink(): raise FileNotFoundError(path)
    expected = {"task_index": index, "pair": pair, "group_index": group_index,
                "candidate_seeds": list(seeds), "hard_low_2d": False,
                "hard_flex": bool(pair_spec["hard_flex"]),
                "hard_scaffold": bool(pair_spec["hard_scaffold"]),
                "preregister_sha256": sha256(PREREG),
                "candidate_pose_sha256": [sha256(path) for path in poses],
                "repair_metadata_sha256": [sha256(path) for path in repair_meta],
                "generation_metadata_sha256": [sha256(path) for path in generation_meta],
                "consensus_weight": CONSENSUS_WEIGHT,
                "mapping_sha256": sha256(pathlib.Path(correspondence.__file__)),
                "worker_sha256": sha256(pathlib.Path(__file__))}
    if audit_only:
        return {**expected, "status": "audit_pass", "selector_reads_native": False,
                "endpoint_metrics_read": False, "threads_max": 1, "gpu_hidden": True}
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or os.getpid() != os.getpgrp():
        raise RuntimeError("GPU or process-group isolation missing")
    task_root = OUTPUT / pair / f"group_{group_index}"
    result_path = task_root / "result.json"
    if result_path.is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if all(payload.get(key) == value for key, value in expected.items()):
            return {**expected, "status": "resumed", "selected_seed": payload["selected_seed"]}
        raise RuntimeError("resume identity drift")
    molecules = [next(m for m in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if m is not None) for path in poses]
    distances, mapping_audit = matrix(molecules)
    consensus = [float(np.median(np.delete(distances[i], i))) for i in range(len(seeds))]
    receptor_scores = []
    pb_pass = []
    for gen_path, rep_path in zip(generation_meta, repair_meta):
        generation = json.loads(gen_path.read_text(encoding="utf-8"))
        repair = json.loads(rep_path.read_text(encoding="utf-8"))
        selected = [item for item in generation["scores"]
                    if item["unique_index"] == generation["selected_unique_index"]]
        if len(selected) != 1: raise RuntimeError("selected generation score mismatch")
        receptor_scores.append(float(selected[0]["selection_score"]))
        pb_pass.append(bool(repair["final_PB_pass"]))
    z_score = zscores(receptor_scores)
    z_consensus = zscores(consensus)
    objective = z_score + CONSENSUS_WEIGHT * z_consensus
    eligible = [i for i, valid in enumerate(pb_pass) if valid]
    if not eligible: eligible = list(range(len(seeds)))
    selected_index = min(eligible, key=lambda i: (float(objective[i]), seeds[i]))
    selected_seed = seeds[selected_index]
    if fault_after_selection: raise SystemExit(97)
    if hold_after_selection > 0: time.sleep(hold_after_selection)
    task_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.",
                                         dir=task_root.parent)); ACTIVE = stage
    shutil.copy2(poses[selected_index], stage / "selected.sdf")
    candidates = [{"seed": seeds[i], "PB20_pass": pb_pass[i],
                   "receptor_score": receptor_scores[i], "consensus_median_A": consensus[i],
                   "z_receptor_score": float(z_score[i]),
                   "z_consensus": float(z_consensus[i]),
                   "objective": float(objective[i])} for i in range(len(seeds))]
    payload = {**expected, "status": "success", "selected_seed": selected_seed,
               "selected_pose_sha256": sha256(stage / "selected.sdf"),
               "candidates": candidates, "chemical_mapping_audit": mapping_audit,
               "distance_definition": "minimum connectivity-preserving full-graph, resonance-aware RMSD; no rigid fit",
               "selector": "PB-valid-first minimum z(receptor_score) + 5*z(consensus_median)",
               "development_scope": "consensus weight fixed using seed group 1 only",
               "selector_reads_native": False, "query_native_ligand_read": False,
               "target_pose_metric_read": False, "SOTA": False}
    atomic_json(stage / "result.json", payload)
    os.replace(stage, task_root); ACTIVE = None
    return {**expected, "status": "computed", "selected_seed": selected_seed}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--fault-after-selection", action="store_true")
    parser.add_argument("--hold-after-selection", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(run(args.index, args.audit_only, args.fault_after_selection,
                         args.hold_after_selection), sort_keys=True))


if __name__ == "__main__": main()
