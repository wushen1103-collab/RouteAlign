#!/usr/bin/env python3
"""Audited PoseBusters evaluation of native-free LS-align group selections."""
from __future__ import annotations
from routealign import settings

import argparse
import atexit
import csv
import hashlib
import json
import os
import pathlib
import shutil
import signal
import tempfile
import time


ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
PREREG = ROOT / "aligndockbench_lsalign_group_pb20_preregister_.md"
GENERATION = ROOT / "runs/lsalign_generate"
ENDPOINT = ROOT / "runs/lsalign_evaluate"
DATA = ROOT / "data/AlignDockBench"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ADB_LSALIGN_GROUP_PB_OUTPUT",
    str(ROOT / "runs/lsalign_validate"),
)).resolve()
EXPECTED_MANIFEST = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
SEEDS = (101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
         151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
         199, 211, 223, 227, 229)
SEED_GROUPS = tuple(SEEDS[offset:offset + 5] for offset in range(0, 25, 5))
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


def cleanup():
    global ACTIVE
    if ACTIVE and ACTIVE.exists():
        shutil.rmtree(ACTIVE)
    ACTIVE = None


def stop(signum, _frame):
    cleanup()
    raise SystemExit(128 + signum)


atexit.register(cleanup)
signal.signal(signal.SIGTERM, stop)
signal.signal(signal.SIGINT, stop)


def cohort():
    if sha256(MANIFEST) != EXPECTED_MANIFEST:
        raise RuntimeError("manifest drift")
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["final_primary_eligible"] == "True"]
    rows.sort(key=lambda row: (row["family"], row["query_id"]))
    if len(rows) != 44:
        raise RuntimeError("cohort drift")
    return rows


def tasks():
    return [(row, group_index, seeds)
            for row in cohort()
            for group_index, seeds in enumerate(SEED_GROUPS, start=1)]


def selected_seed(pair, seeds):
    candidates = []
    hashes = []
    for seed in seeds:
        path = GENERATION / pair / f"seed_{seed}" / "metadata.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        if (payload.get("status") != "success" or payload.get("pair") != pair or
                payload.get("seed") != seed or payload.get("query_native_ligand_read") is not False or
                payload.get("selector_reads_native") is not False):
            raise RuntimeError("generation identity mismatch: " + str(path))
        metric = payload["selected_metrics"]
        key = (-min(metric["pc_score_query"], metric["pc_score_template"]),
               -metric["jaccard"], metric["ls_rmsd_A"], metric["mode"], seed)
        candidates.append((key, seed, path))
        hashes.append(sha256(path))
    _key, seed, path = min(candidates)
    return seed, path, hashes


def run(index, fault_after_workload=False, hold_after_workload=0.0, audit_only=False):
    global ACTIVE
    if pathlib.Path.cwd().resolve() != ROOT.resolve():
        raise RuntimeError("wrong project directory")
    all_tasks = tasks()
    if not 0 <= index < len(all_tasks) or (os.cpu_count() or 0) < 1:
        raise RuntimeError("scope/resources")
    row, group_index, seeds = all_tasks[index]
    pair = f"{row['family']}__{row['query_id']}"
    seed, metadata_path, group_metadata_hashes = selected_seed(pair, seeds)
    pose = GENERATION / pair / f"seed_{seed}" / "selected.sdf"
    protein = DATA / row["family"] / row["template_id"] / "protein.pdb"
    endpoint_path = ENDPOINT / pair / f"seed_{seed}" / "result.json"
    for path in (pose, protein, endpoint_path, metadata_path):
        path.resolve().relative_to(ROOT.resolve())
        if not path.is_file():
            raise FileNotFoundError(path)
    endpoint = json.loads(endpoint_path.read_text(encoding="utf-8"))
    if (endpoint.get("status") != "success" or endpoint.get("pair") != pair or
            endpoint.get("seed") != seed or endpoint.get("generator_reads_native") is not False):
        raise RuntimeError("endpoint identity mismatch")
    expected = {
        "task_index": index,
        "pair": pair,
        "group_index": group_index,
        "group_seeds": list(seeds),
        "selected_seed": seed,
        "selection_rule": "max min(PC-scoreQ,PC-scoreT), then max Jaccard, then min LS-RMSD, then fixed mode and seed",
        "selector_reads_native": False,
        "manifest_sha256": EXPECTED_MANIFEST,
        "preregister_sha256": sha256(PREREG),
        "pose_sha256": sha256(pose),
        "protein_sha256": sha256(protein),
        "endpoint_sha256": sha256(endpoint_path),
        "generation_metadata_sha256": sha256(metadata_path),
        "group_generation_metadata_sha256": group_metadata_hashes,
        "worker_sha256": sha256(__file__),
    }
    if audit_only:
        return {**expected, "status": "audit_pass", "pose_bytes": pose.stat().st_size,
                "protein_bytes": protein.stat().st_size}
    task = OUTPUT / pair / f"group_{group_index}"
    result_path = task / "result.json"
    if result_path.is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if all(payload.get(key) == value for key, value in expected.items()):
            return {**expected, "status": "resumed", "PB20_pass": payload["PB20_pass"]}
        raise RuntimeError("resume drift")
    from posebusters import PoseBusters, __version__
    task.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task.name + ".staging.", dir=task.parent))
    ACTIVE = stage
    frame = PoseBusters(config="dock", max_workers=0).bust(
        str(pose), mol_cond=str(protein), full_report=False)
    if len(frame) != 1:
        raise RuntimeError("PoseBusters row count")
    checks = {str(key): bool(value) for key, value in frame.iloc[0].to_dict().items()}
    scientific = [key for key in checks if key not in {"mol_pred_loaded", "mol_cond_loaded"}]
    if len(scientific) != 19:
        raise RuntimeError(f"expected 19 scientific checks, got {len(scientific)}")
    if fault_after_workload:
        cleanup()
        raise SystemExit(97)
    if hold_after_workload > 0:
        time.sleep(hold_after_workload)
    payload = {
        **expected,
        "status": "success",
        "posebusters_version": __version__,
        "config": "dock",
        "target_native_ligand_read": False,
        "template_receptor_used_for_physical_checks": True,
        "rmsd_A": float(endpoint["rmsd_A"]),
        "success_le_2A": bool(endpoint["success_le_2A"]),
        "PB_selected_checks": len(scientific),
        "PB_pass_count": sum(checks[key] for key in scientific),
        "PB20_pass": all(checks.values()),
        "PB_scientific_pass": all(checks[key] for key in scientific),
        "joint_success": bool(endpoint["success_le_2A"] and all(checks.values())),
        "failed": [key for key in scientific if not checks[key]],
        "checks": checks,
        "SOTA": False,
    }
    atomic_json(stage / "result.json", payload)
    os.replace(stage, task)
    ACTIVE = None
    return {**expected, "status": "computed", "PB20_pass": payload["PB20_pass"],
            "joint_success": payload["joint_success"], "failed": payload["failed"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--fault-after-workload", action="store_true")
    parser.add_argument("--hold-after-workload", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(run(args.index, args.fault_after_workload,
                         args.hold_after_workload, args.audit_only), sort_keys=True))


if __name__ == "__main__":
    main()
