#!/usr/bin/env python3
"""Controlled same-row ROSHAMBO (JCIM 2024) generator for AlignDockBench full44.

The generator sees the public template ligand and a coordinate-free query SMILES.
It never opens the native query pose.  One task is one pair/seed and commits through
an atomic directory rename.  The ROSHAMBO multiprocessing helpers are capped at the
current eight-CPU affinity mask before any pool is created.
"""
from __future__ import annotations
from routealign import settings

import argparse
import atexit
import csv
import fcntl
import hashlib
import json
import os
import pathlib
import shutil
import signal
import subprocess
import tempfile
import time

from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import rdMolAlign, rdMolTransforms


ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
PREREG = ROOT / "config/protocol.json"
DATASET = ROOT / "data/AlignDockBench"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ROSHAMBO_OUTPUT",
    str(ROOT / "runs/roshambo_generate"),
)).resolve()
SOURCE = ROOT / "third_party/roshambo"
EXPECTED_MANIFEST_SHA256 = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
EXPECTED_SOURCE_COMMIT = "916fe97f81a3fdb48e7bb85384eb597836d2b01e"
EXPECTED_SETUP_SHA256 = "d4f1ae0eea155a790db17f174a5784db77e60d10e73c07fbf1a2efa121936a7f"
EXPECTED_EXTENSION_SHA256 = "9f4d8899452e89abc1b441fd1453b16e5889ce4548bbd2de40a22e4f90308830"
SEEDS = (101, 127, 151, 179, 199)
PHYSICAL_GPU = int(os.environ.get("ROUTEALIGN_GPU", "-1"))
GPU_MAX_MEMORY_MIB = 500
GPU_MAX_UTIL_PERCENT = 10
CPU_CAP = 8
ACTIVE: pathlib.Path | None = None
LOCK_HANDLE = None
OWNER_PID = os.getpid()


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
    global ACTIVE, LOCK_HANDLE
    # multiprocessing uses fork here.  Pool children inherit atexit handlers and
    # lock descriptors, but must never clean the parent task's staging directory.
    if os.getpid() != OWNER_PID:
        return
    if ACTIVE is not None and ACTIVE.exists():
        shutil.rmtree(ACTIVE)
    ACTIVE = None
    if LOCK_HANDLE is not None:
        try:
            fcntl.flock(LOCK_HANDLE, fcntl.LOCK_UN)
            LOCK_HANDLE.close()
        except Exception:
            pass
        LOCK_HANDLE = None


def stop(sig, _frame) -> None:
    if os.getpid() != OWNER_PID:
        raise SystemExit(128 + sig)
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
        rows = [r for r in csv.DictReader(handle) if r["final_primary_eligible"] == "True"]
    by_pair = {f"{r['family']}__{r['query_id']}": r for r in rows}
    pairs = list(prereg["pairs"])
    if len(pairs) != 44 or any(p["pair"] not in by_pair for p in pairs):
        raise RuntimeError("full44 cohort drift")
    return [(i, item, by_pair[item["pair"]], seed) for i, (item, seed) in enumerate(
        (item, seed) for item in pairs for seed in SEEDS
    )]


def source_identity() -> dict:
    commit = subprocess.check_output(
        ["git", "-C", str(SOURCE), "rev-parse", "HEAD"], text=True
    ).strip()
    setup_sha = sha256(SOURCE / "setup.py")
    extension = next((SOURCE / "roshambo").glob("cpaper*.so"))
    ext_sha = sha256(extension)
    if commit != EXPECTED_SOURCE_COMMIT or setup_sha != EXPECTED_SETUP_SHA256:
        raise RuntimeError("ROSHAMBO source/build identity drift")
    return {"roshambo_git_commit": commit, "setup_py_sha256": setup_sha,
            "cpaper_extension_sha256": ext_sha}


def gpu_state() -> dict:
    line = subprocess.check_output([
        "nvidia-smi", f"--id={PHYSICAL_GPU}",
        "--query-gpu=index,uuid,memory.used,memory.total,utilization.gpu,pstate",
        "--format=csv,noheader,nounits",
    ], text=True, timeout=15).strip()
    fields = [x.strip() for x in line.split(",")]
    if len(fields) != 6 or int(fields[0]) != PHYSICAL_GPU:
        raise RuntimeError("cannot identify physical GPU")
    state = {"index": int(fields[0]), "uuid": fields[1], "memory_used_mib": int(fields[2]),
             "memory_total_mib": int(fields[3]), "utilization_percent": int(fields[4]),
             "pstate": fields[5]}
    if state["memory_used_mib"] >= GPU_MAX_MEMORY_MIB or state["utilization_percent"] > GPU_MAX_UTIL_PERCENT:
        raise RuntimeError("requested GPU no longer satisfies preregistered idle gate: " + json.dumps(state))
    return state


def acquire_gpu_lock() -> None:
    global LOCK_HANDLE
    lock_dir = ROOT / ".locks"
    lock_dir.mkdir(exist_ok=True)
    LOCK_HANDLE = (lock_dir / f"gpu{PHYSICAL_GPU}_roshambo.lock").open("w")
    try:
        fcntl.flock(LOCK_HANDLE, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise RuntimeError("ROSHAMBO device lock already held") from exc
    LOCK_HANDLE.write(f"pid={os.getpid()}\n")
    LOCK_HANDLE.flush()


def spec_for(index: int) -> tuple[dict, dict, pathlib.Path]:
    tasks = table()
    if not 0 <= index < len(tasks):
        raise IndexError(index)
    _, item, row, seed = tasks[index]
    template = DATASET / row["family"] / row["template_id"] / "ligand.sdf"
    template.resolve().relative_to(ROOT)
    if not template.is_file() or template.is_symlink():
        raise FileNotFoundError(template)
    identity = source_identity()
    spec = {
        "task_index": index, "pair": item["pair"], "family": row["family"],
        "template_id": row["template_id"], "query_id": row["query_id"], "seed": seed,
        "template_sha256": sha256(template), "query_graph_sha256": row["query_graph_sha256"],
        "manifest_sha256": EXPECTED_MANIFEST_SHA256, "preregister_sha256": sha256(PREREG),
        "generator_reads_native_query_pose": False, "query_source": "manifest query_input_smiles",
        "reference_source": "public template ligand.sdf", "n_confs_requested": 10,
        "conformer_method": "ETKDGv2", "conformer_force_field": "MMFF94s",
        "optimize_conformers": True, "ignore_hydrogens_during_alignment": True,
        "volume_type": "gaussian", "color": True, "sort_by": "ComboTanimoto",
        "max_conformers_written": 1, "physical_gpu": PHYSICAL_GPU,
        "visible_gpu": 0, "cpu_affinity_required": CPU_CAP, "rdkit_version": rdBase.rdkitVersion,
        **identity,
    }
    return spec, row, template


def valid(result_dir: pathlib.Path, spec: dict) -> bool:
    metadata = result_dir / "metadata.json"
    selected = result_dir / "selected.sdf"
    if not metadata.is_file() or not selected.is_file():
        return False
    try:
        payload = json.loads(metadata.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return (payload.get("status") == "success"
            and all(payload.get(k) == v for k, v in spec.items())
            and payload.get("selected_sha256") == sha256(selected))


def transform_back(hits: pathlib.Path, template_path: pathlib.Path, selected_path: pathlib.Path) -> tuple[int, float]:
    hits_mols = [m for m in Chem.SDMolSupplier(str(hits), removeHs=False, sanitize=True) if m is not None]
    if len(hits_mols) != 2:
        raise RuntimeError(f"expected reference plus top-1 hit, found {len(hits_mols)}")
    output_ref, selected = hits_mols
    template = next((m for m in Chem.SDMolSupplier(str(template_path), removeHs=False, sanitize=True) if m is not None), None)
    if template is None:
        raise RuntimeError("invalid template SDF")
    ref_heavy = Chem.RemoveHs(output_ref)
    template_heavy = Chem.RemoveHs(template)
    if (ref_heavy.GetNumAtoms() != template_heavy.GetNumAtoms()
            or [a.GetAtomicNum() for a in ref_heavy.GetAtoms()] != [a.GetAtomicNum() for a in template_heavy.GetAtoms()]):
        raise RuntimeError("reference atom order drift prevents frame restoration")
    rms, matrix = rdMolAlign.GetAlignmentTransform(
        ref_heavy, template_heavy,
        atomMap=[(i, i) for i in range(ref_heavy.GetNumAtoms())],
    )
    rdMolTransforms.TransformConformer(selected.GetConformer(), matrix)
    selected = Chem.RemoveHs(selected)
    writer = Chem.SDWriter(str(selected_path))
    writer.write(selected)
    writer.close()
    return len(hits_mols) - 1, float(rms)


def run(index: int, audit_only: bool = False, fault_after_workload: bool = False,
        hold_after_workload: float = 0.0) -> dict:
    global ACTIVE
    if pathlib.Path.cwd().resolve() != ROOT:
        raise RuntimeError("fixed project root required")
    spec, row, template = spec_for(index)
    if audit_only:
        return {**spec, "status": "audit_pass", "task_count": len(table())}
    if os.environ.get("CUDA_VISIBLE_DEVICES") != str(PHYSICAL_GPU):
        raise RuntimeError("CUDA_VISIBLE_DEVICES must select requested GPU")
    if os.getpid() != os.getpgrp():
        raise RuntimeError("dedicated process group required")
    affinity = sorted(os.sched_getaffinity(0))
    if len(affinity) > CPU_CAP:
        raise RuntimeError(f"CPU affinity exceeds {CPU_CAP}: {affinity}")
    result_dir = OUTPUT / spec["pair"] / f"seed_{spec['seed']}"
    if valid(result_dir, spec):
        return {**json.loads((result_dir / "metadata.json").read_text()), "status": "resumed"}
    result_dir.parent.mkdir(parents=True, exist_ok=True)
    acquire_gpu_lock()
    preflight = gpu_state()
    stage = pathlib.Path(tempfile.mkdtemp(prefix=result_dir.name + ".staging.", dir=result_dir.parent))
    ACTIVE = stage
    work = stage / "work"
    work.mkdir()
    shutil.copyfile(template, work / "template.sdf")
    (work / "query.smi").write_text(f"SMILES Name\n{row['query_input_smiles']} query\n", encoding="utf-8")

    # Import only after device, lock, and affinity checks.  Patch imported aliases
    # because upstream otherwise creates an unrestricted worker pool in three pools.
    import roshambo.core as roshambo_core
    import roshambo.utilities as roshambo_utilities
    from roshambo.api import get_similarity_scores
    pool_size = min(CPU_CAP, len(os.sched_getaffinity(0)))
    roshambo_core.cpu_count = lambda: pool_size
    roshambo_utilities.cpu_count = lambda: pool_size
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib-roshambo"))
    os.environ.setdefault("RDKIT_DATA_DIR", RDConfig.RDDataDir)

    started = time.time()
    get_similarity_scores(
        ref_file="template.sdf", dataset_files_pattern="query.smi",
        ignore_hs=True, n_confs=10, keep_mol=False, working_dir=str(work),
        name_prefix="query", smiles_kwargs={"delimiter": " ", "titleLine": True},
        embed_kwargs={"randomSeed": spec["seed"]}, gpu_id=0,
        volume_type="gaussian", use_carbon_radii=True, color=True,
        sort_by="ComboTanimoto", write_to_file=True, max_conformers=1,
        filename="hits.sdf", random_seed=spec["seed"], method="ETKDGv2",
        ff="MMFF94s", add_hs=True, opt_confs=True, calc_energy=True,
        energy_iters=200, num_threads=1,
    )
    runtime = time.time() - started
    generated = sum(1 for m in Chem.SDMolSupplier(str(work / "mols.sdf"), removeHs=False) if m is not None)
    _, frame_rms = transform_back(work / "hits.sdf", template, stage / "selected.sdf")
    shutil.copyfile(work / "roshambo.csv", stage / "roshambo.tsv")
    if fault_after_workload:
        raise SystemExit(97)
    if hold_after_workload > 0:
        time.sleep(hold_after_workload)
    metadata = {
        **spec, "status": "success", "gpu_preflight": preflight,
        "cpu_affinity": affinity, "pool_size": pool_size,
        "generated_conformers_observed": generated,
        "reference_frame_restore_rmsd_A": frame_rms,
        "runtime_seconds": runtime,
        "selected_sha256": sha256(stage / "selected.sdf"),
        "score_table_sha256": sha256(stage / "roshambo.tsv"),
    }
    atomic_json(stage / "metadata.json", metadata)
    if result_dir.exists():
        if valid(result_dir, spec):
            shutil.rmtree(stage)
            ACTIVE = None
            return {**json.loads((result_dir / "metadata.json").read_text()), "status": "resumed_after_race"}
        raise RuntimeError("refusing to overwrite invalid existing task directory")
    os.replace(stage, result_dir)
    ACTIVE = None
    clean()
    return metadata


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
