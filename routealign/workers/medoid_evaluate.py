#!/usr/bin/env python3
"""Audited endpoint evaluation for the frozen full 44-pair expansion."""
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

from rdkit import Chem, rdBase
from posebusters import __version__ as posebusters_version
from posebusters.modules.rmsd import check_rmsd


ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
DATASET = ROOT / "data/AlignDockBench"
SELECTOR = ROOT / "runs/medoid_select"
REPAIR = ROOT / "runs/bounded_repair"
PREREG = ROOT / "config/protocol.json"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ADB_MEDOID_EVAL_OUTPUT",
    str(ROOT / "runs/medoid_evaluate"),
)).resolve()
EXPECTED_MANIFEST_SHA256 = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
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


def atomic_json(path, payload):
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def task_table():
    if sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError("manifest identity drift")
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    pair_specs = prereg["pairs"]
    pairs = [item["pair"] for item in pair_specs]
    groups = prereg["outer_seed_groups"]
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["final_primary_eligible"] == "True"]
    by_pair = {f"{row['family']}__{row['query_id']}": row for row in rows}
    if set(pairs) - set(by_pair):
        raise RuntimeError("pair/manifest mismatch")
    return [(index, pair_spec, group_index, tuple(group), by_pair[pair_spec["pair"]])
            for index, (pair_spec, group_index, group) in enumerate(
                (p, gi, g) for p in pair_specs for gi, g in enumerate(groups, start=1)
            )]


def load_single(path):
    molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if mol is not None]
    if len(molecules) != 1 or molecules[0].GetNumConformers() != 1:
        raise RuntimeError("pose parse failed: " + str(path))
    return molecules[0]


def run(index, audit_only=False, fault_after_metric=False, hold_after_metric=0.0):
    global ACTIVE
    table = task_table()
    if pathlib.Path.cwd().resolve() != ROOT.resolve() or not 0 <= index < len(table):
        raise RuntimeError("scope/index mismatch")
    _, pair_spec, group_index, group, row = table[index]
    pair = pair_spec["pair"]
    pose = SELECTOR / pair / f"group_{group_index}" / "selected.sdf"
    selection_path = SELECTOR / pair / f"group_{group_index}" / "result.json"
    native = DATASET / row["family"] / row["query_id"] / "ligand.sdf"
    for path in (pose, selection_path, native):
        path.resolve().relative_to(ROOT.resolve())
        if not path.is_file():
            raise FileNotFoundError(path)
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    selected_seed = int(selection["selected_seed"])
    if selected_seed not in group:
        raise RuntimeError("selected seed outside frozen group")
    repair_path = REPAIR / pair / f"seed_{selected_seed}" / "result.json"
    if not repair_path.is_file():
        raise FileNotFoundError(repair_path)
    repair = json.loads(repair_path.read_text(encoding="utf-8"))
    expected = {
        "task_index": index,
        "pair": pair,
        "group_index": group_index,
        "selected_seed": selected_seed,
        "hard_low_2d": bool(pair_spec["hard_low_2d"]),
        "hard_flex": bool(pair_spec["hard_flex"]),
        "hard_scaffold": bool(pair_spec["hard_scaffold"]),
        "manifest_sha256": EXPECTED_MANIFEST_SHA256,
        "preregister_sha256": sha256(PREREG),
        "prediction_sha256": sha256(pose),
        "selector_result_sha256": sha256(selection_path),
        "repair_result_sha256": sha256(repair_path),
        "native_sha256": sha256(native),
        "worker_sha256": sha256(__file__),
    }
    if audit_only:
        return {**expected, "status": "audit_pass", "endpoint_only": True,
                "threads_max": 1, "gpu_hidden": True}
    task_root = OUTPUT / pair / f"group_{group_index}"
    result_path = task_root / "result.json"
    if result_path.is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if all(payload.get(key) == value for key, value in expected.items()):
            return {**expected, "status": "resumed", "rmsd_A": payload["rmsd_A"]}
        raise RuntimeError("resume identity drift")
    metric = check_rmsd(load_single(pose), load_single(native), rmsd_threshold=2.0,
                        heavy_only=True, choose_by="rmsd")["results"]
    rmsd = float(metric["rmsd"])
    if fault_after_metric:
        raise SystemExit(97)
    if hold_after_metric > 0:
        time.sleep(hold_after_metric)
    task_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.", dir=task_root.parent))
    ACTIVE = stage
    success = bool(metric["rmsd_within_threshold"])
    pb_pass = bool(repair["PB20_pass"])
    payload = {
        **expected,
        "status": "success",
        "rmsd_A": rmsd,
        "success_le_2A": success,
        "PB20_pass": pb_pass,
        "joint_success": bool(success and pb_pass),
        "repair_decision": repair["previous_decision"],
        "trust125_override": bool(repair["override"]),
        "query_native_ligand_read_for_endpoint": True,
        "selector_read_native": False,
        "posebusters_version": posebusters_version,
        "rdkit_version": rdBase.rdkitVersion,
        "SOTA": False,
    }
    atomic_json(stage / "result.json", payload)
    os.replace(stage, task_root)
    ACTIVE = None
    return {**expected, "status": "computed", "rmsd_A": rmsd,
            "success_le_2A": success, "PB20_pass": pb_pass,
            "joint_success": bool(success and pb_pass)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--fault-after-metric", action="store_true")
    parser.add_argument("--hold-after-metric", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(run(args.index, args.audit_only, args.fault_after_metric,
                         args.hold_after_metric), sort_keys=True))


if __name__ == "__main__":
    main()
