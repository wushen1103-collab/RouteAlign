#!/usr/bin/env python3
"""Single-process PoseBusters evaluation for frozen ROSHAMBO full44 outputs."""
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
PREREG = ROOT / "config/roshambo_protocol.md"
FULL44_PREREG = ROOT / "config/protocol.json"
GENERATION = ROOT / "runs/roshambo_generate"
ENDPOINT = ROOT / "runs/roshambo_evaluate"
DATA = ROOT / "data/AlignDockBench"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ROSHAMBO_PB20_OUTPUT",
    str(ROOT / "runs/roshambo_validate"),
)).resolve()
EXPECTED_MANIFEST = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
SEEDS = (101, 127, 151, 179, 199)
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


def tasks():
    if sha256(MANIFEST) != EXPECTED_MANIFEST:
        raise RuntimeError("manifest drift")
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["final_primary_eligible"] == "True"]
    by_pair = {f"{r['family']}__{r['query_id']}": r for r in rows}
    pairs = json.loads(FULL44_PREREG.read_text(encoding="utf-8"))["pairs"]
    if len(pairs) != 44 or any(p["pair"] not in by_pair for p in pairs):
        raise RuntimeError("cohort drift")
    return [(item, by_pair[item["pair"]], seed) for item in pairs for seed in SEEDS]


def run(index, fault_after_workload=False, hold_after_workload=0.0, audit_only=False):
    global ACTIVE
    if pathlib.Path.cwd().resolve() != ROOT:
        raise RuntimeError("wrong project directory")
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or os.getpid() != os.getpgrp():
        raise RuntimeError("GPU/process-group isolation missing")
    if len(os.sched_getaffinity(0)) != 1:
        raise RuntimeError("PoseBusters must run on exactly one affine CPU")
    all_tasks = tasks()
    if not 0 <= index < len(all_tasks):
        raise IndexError(index)
    item, row, seed = all_tasks[index]
    pair = item["pair"]
    pose = GENERATION / pair / f"seed_{seed}" / "selected.sdf"
    generation_metadata = GENERATION / pair / f"seed_{seed}" / "metadata.json"
    protein = DATA / row["family"] / row["template_id"] / "protein.pdb"
    endpoint_path = ENDPOINT / pair / f"seed_{seed}" / "result.json"
    for path in (pose, generation_metadata, protein, endpoint_path):
        path.resolve().relative_to(ROOT)
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(path)
    generated = json.loads(generation_metadata.read_text(encoding="utf-8"))
    endpoint = json.loads(endpoint_path.read_text(encoding="utf-8"))
    if (generated.get("status") != "success" or generated.get("pair") != pair
            or generated.get("seed") != seed
            or generated.get("generator_reads_native_query_pose") is not False):
        raise RuntimeError("generation identity mismatch")
    if (endpoint.get("status") != "success" or endpoint.get("pair") != pair
            or endpoint.get("seed") != seed or endpoint.get("generator_reads_native") is not False):
        raise RuntimeError("endpoint identity mismatch")
    expected = {
        "task_index": index, "pair": pair, "seed": seed,
        "manifest_sha256": EXPECTED_MANIFEST,
        "roshambo_preregister_sha256": sha256(PREREG),
        "full44_preregister_sha256": sha256(FULL44_PREREG),
        "pose_sha256": sha256(pose), "protein_sha256": sha256(protein),
        "generation_metadata_sha256": sha256(generation_metadata),
        "endpoint_sha256": sha256(endpoint_path), "generator_reads_native": False,
        "target_native_ligand_read_for_physical_checks": False,
        "template_receptor_used_for_physical_checks": True,
        "worker_sha256": sha256(__file__),
    }
    if audit_only:
        return {**expected, "status": "audit_pass", "task_count": len(all_tasks)}
    task = OUTPUT / pair / f"seed_{seed}"
    result_path = task / "result.json"
    if result_path.is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if all(payload.get(key) == value for key, value in expected.items()):
            return {**expected, "status": "resumed", "PB20_pass": payload["PB20_pass"],
                    "joint_success": payload["joint_success"]}
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
        **expected, "status": "success", "posebusters_version": __version__,
        "config": "dock", "rmsd_A": float(endpoint["rmsd_A"]),
        "success_le_2A": bool(endpoint["success_le_2A"]),
        "PB_selected_checks": len(scientific),
        "PB_pass_count": sum(checks[key] for key in scientific),
        "PB20_pass": all(checks.values()),
        "PB_scientific_pass": all(checks[key] for key in scientific),
        "joint_success": bool(endpoint["success_le_2A"] and all(checks.values())),
        "failed": [key for key in scientific if not checks[key]],
        "checks": checks, "SOTA": False,
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
