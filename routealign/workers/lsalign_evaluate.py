#!/usr/bin/env python3
"""Endpoint-only evaluation for frozen full44 LS-align predictions."""
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
PREREG = ROOT / "config/protocol.json"
DATASET = ROOT / "data/AlignDockBench"
PREDICTIONS = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ENDPOINT_PREDICTIONS",
    str(ROOT / "runs/lsalign_generate"),
)).resolve()
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ENDPOINT_OUTPUT",
    os.environ.get("ROUTEALIGN_LSALIGN_ENDPOINT_OUTPUT",
                   str(ROOT / "runs/lsalign_evaluate")),
)).resolve()
EXPECTED_MANIFEST_SHA256 = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
SEEDS = (101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163,
         167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229)
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


def table() -> list[tuple[int, dict, dict, int]]:
    if sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError("manifest identity drift")
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["final_primary_eligible"] == "True"]
    by_pair = {f"{row['family']}__{row['query_id']}": row for row in rows}
    pairs = list(prereg["pairs"])
    if len(pairs) != 44 or any(item["pair"] not in by_pair for item in pairs):
        raise RuntimeError("full44 cohort drift")
    return [(index, item, by_pair[item["pair"]], seed) for index, (item, seed) in enumerate(
        (item, seed) for item in pairs for seed in SEEDS
    )]


def load_single(path: pathlib.Path) -> Chem.Mol:
    molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if mol is not None]
    if len(molecules) != 1 or molecules[0].GetNumConformers() != 1:
        raise RuntimeError("invalid SDF: " + str(path))
    return molecules[0]


def expected(index: int) -> tuple[dict, pathlib.Path, pathlib.Path]:
    tasks = table()
    if not 0 <= index < len(tasks):
        raise IndexError(index)
    _, item, row, seed = tasks[index]
    prediction = PREDICTIONS / item["pair"] / f"seed_{seed}" / "selected.sdf"
    metadata = PREDICTIONS / item["pair"] / f"seed_{seed}" / "metadata.json"
    native = DATASET / row["family"] / row["query_id"] / "ligand.sdf"
    for path in (prediction, metadata, native):
        path.resolve().relative_to(ROOT)
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(path)
    spec = {"task_index": index, "pair": item["pair"], "seed": seed,
            "prediction_root": str(PREDICTIONS),
            "prediction_sha256": sha256(prediction), "prediction_metadata_sha256": sha256(metadata),
            "native_sha256": sha256(native), "manifest_sha256": EXPECTED_MANIFEST_SHA256,
            "preregister_sha256": sha256(PREREG), "endpoint_reads_native": True,
            "generator_reads_native": False, "threads_max": 1, "gpu_hidden": True,
            "rdkit_version": rdBase.rdkitVersion, "posebusters_version": posebusters_version}
    return spec, prediction, native


def valid(path: pathlib.Path, spec: dict) -> bool:
    if not path.is_file():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return payload.get("status") == "success" and all(payload.get(k) == v for k, v in spec.items())


def run(index: int, audit_only: bool = False, fault_after_metric: bool = False,
        hold_after_metric: float = 0.0) -> dict:
    global ACTIVE
    if pathlib.Path.cwd().resolve() != ROOT:
        raise RuntimeError("fixed project root required")
    spec, prediction, native = expected(index)
    if audit_only:
        return {**spec, "status": "audit_pass"}
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or os.getpid() != os.getpgrp():
        raise RuntimeError("GPU or process-group isolation missing")
    result_root = OUTPUT / spec["pair"] / f"seed_{spec['seed']}"
    result_path = result_root / "result.json"
    if valid(result_path, spec):
        return {**json.loads(result_path.read_text()), "status": "resumed"}
    metric = check_rmsd(load_single(prediction), load_single(native), rmsd_threshold=2.0,
                        heavy_only=True, choose_by="rmsd")["results"]
    rmsd = float(metric["rmsd"])
    if fault_after_metric:
        raise SystemExit(97)
    if hold_after_metric > 0:
        time.sleep(hold_after_metric)
    result_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=result_root.name + ".staging.", dir=result_root.parent))
    ACTIVE = stage
    payload = {**spec, "status": "success", "rmsd_A": rmsd,
               "success_le_2A": bool(metric["rmsd_within_threshold"]), "SOTA": False}
    atomic_json(stage / "result.json", payload)
    os.replace(stage, result_root)
    ACTIVE = None
    return payload


def main() -> None:
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
