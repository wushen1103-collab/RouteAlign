#!/usr/bin/env python3
"""Materialize the outcome-blind 1.25 A repair policy for all 44 pairs."""
from __future__ import annotations
from routealign import settings

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


ROOT = settings.workspace()
PREREG = ROOT / "config/protocol.json"
INPUT = ROOT / "runs/repair"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ADB_TRUST125_OUTPUT",
    str(ROOT / "runs/bounded_repair"),
)).resolve()
ACTIVE = None


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


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


def table():
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    seeds = [seed for group in prereg["outer_seed_groups"] for seed in group]
    if len(prereg["pairs"]) != 44 or len(seeds) != 25:
        raise RuntimeError("preregistration mismatch")
    return [(index, item, seed) for index, (item, seed) in enumerate(
        (item, seed) for item in prereg["pairs"] for seed in seeds)]


def run(index, audit_only=False, fault_after_decision=False, hold_after_decision=0.0):
    global ACTIVE
    tasks = table()
    if pathlib.Path.cwd().resolve() != ROOT.resolve() or not 0 <= index < len(tasks):
        raise RuntimeError("scope/index mismatch")
    _, item, seed = tasks[index]
    pair = item["pair"]
    base = INPUT / pair / f"seed_{seed}"
    old_result_path = base / "result.json"
    old_selected = base / "selected.sdf"
    for path in (old_result_path, old_selected):
        if not path.is_file():
            raise FileNotFoundError(path)
    old = json.loads(old_result_path.read_text(encoding="utf-8"))
    candidate = base / "candidate.sdf"
    override = (old["decision"] == "repair_rejected" and float(old["movement_A"]) <= 1.25
                and not old["repaired_failed"] and candidate.is_file())
    source = candidate if override else old_selected
    pb_pass = True if override else bool(old["final_PB_pass"])
    expected = {"task_index": index, "pair": pair, "seed": seed,
                "preregister_sha256": sha256(PREREG), "old_result_sha256": sha256(old_result_path),
                "source_pose_sha256": sha256(source), "source_role": "candidate" if override else "selected",
                "worker_sha256": sha256(__file__)}
    if audit_only:
        return {**expected, "status": "audit_pass", "override": override,
                "threads_max": 1, "gpu_hidden": True}
    task_root = OUTPUT / pair / f"seed_{seed}"
    result_path = task_root / "result.json"
    if result_path.is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if all(payload.get(key) == value for key, value in expected.items()):
            return {**expected, "status": "resumed", "override": payload["override"]}
        raise RuntimeError("resume identity drift")
    if fault_after_decision:
        raise SystemExit(97)
    if hold_after_decision > 0:
        time.sleep(hold_after_decision)
    task_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.", dir=task_root.parent))
    ACTIVE = stage
    shutil.copy2(source, stage / "selected.sdf")
    payload = {**expected, "status": "success", "override": override,
               "trust_region_A": 1.25, "PB20_pass": pb_pass,
               "movement_A": float(old["movement_A"]), "previous_decision": old["decision"],
               "selected_sha256": sha256(stage / "selected.sdf"),
               "query_native_ligand_read": False, "target_pose_metric_read": False, "SOTA": False}
    (stage / "result.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(stage, task_root)
    ACTIVE = None
    return {**expected, "status": "computed", "override": override, "PB20_pass": pb_pass}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--fault-after-decision", action="store_true")
    parser.add_argument("--hold-after-decision", type=float, default=0.0)
    args = parser.parse_args()
    if args.list:
        for index, item, seed in table():
            print(index, item["pair"], seed)
        return
    if args.index is None:
        parser.error("--index required")
    print(json.dumps(run(args.index, args.audit_only, args.fault_after_decision,
                         args.hold_after_decision), sort_keys=True))


if __name__ == "__main__":
    main()
