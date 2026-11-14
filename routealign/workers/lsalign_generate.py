#!/usr/bin/env python3
"""Outcome-blind LS-align expert on all 44 preregistered conditional pairs."""
from __future__ import annotations
from routealign import settings

import argparse
import atexit
import csv
import hashlib
import json
import math
import os
import pathlib
import resource
import shutil
import signal
import subprocess
import tempfile
import time

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D


ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
PREREG = ROOT / "config/protocol.json"
DATASET = ROOT / "data/AlignDockBench"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_LSALIGN_EXPERT_OUTPUT",
    str(ROOT / "runs/lsalign_generate"),
)).resolve()
LSALIGN = ROOT / "tools/LSalign"
OBABEL = settings.openbabel()
EXPECTED_MANIFEST_SHA256 = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
EXPECTED_LSALIGN_SHA256 = "598afc0cd5f8994d9aba159b45315de12fa5ecdff4a226de07b281959c9bf1bc"
EXPECTED_OBABEL_SHA256 = "22d18273dc9277ac174be21eb5bfd201ca9c364d03f00d1c47d7960da46032a3"
SEEDS = (101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163,
         167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229)
MODES = (("rigid_acc", 0), ("flex_acc", 1))
CPU = settings.cpu_index()
ACTIVE: pathlib.Path | None = None


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def scope(path: pathlib.Path, must_exist: bool = True) -> pathlib.Path:
    resolved = path.resolve()
    resolved.relative_to(ROOT)
    if path.is_symlink():
        raise RuntimeError("symlink forbidden in controlled project path: " + str(path))
    if must_exist and not path.exists():
        raise FileNotFoundError(path)
    return resolved


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


def atomic_json(path: pathlib.Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def load_rows() -> tuple[list[dict], dict]:
    if sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError("manifest identity drift")
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["final_primary_eligible"] == "True"]
    by_pair = {f"{row['family']}__{row['query_id']}": row for row in rows}
    specs = list(prereg["pairs"])
    if len(specs) != 44 or any(item["pair"] not in by_pair for item in specs):
        raise RuntimeError("full44 cohort drift")
    return [{**item, "row": by_pair[item["pair"]]} for item in specs], prereg


def tasks() -> list[tuple[int, dict, int]]:
    specs, _ = load_rows()
    return [(index, item, seed) for index, (item, seed) in enumerate(
        (item, seed) for item in specs for seed in SEEDS
    )]


def load_single(path: pathlib.Path) -> Chem.Mol:
    molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if mol is not None]
    if len(molecules) != 1 or molecules[0].GetNumConformers() != 1:
        raise RuntimeError("invalid single-conformer SDF: " + str(path))
    return molecules[0]


def graph_identity(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def embed_query(smiles: str, seed: int) -> tuple[Chem.Mol, dict]:
    base = Chem.MolFromSmiles(smiles)
    if base is None:
        raise RuntimeError("query SMILES parse failed")
    molecule = Chem.AddHs(base)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.numThreads = 1
    params.useRandomCoords = False
    code = int(AllChem.EmbedMolecule(molecule, params))
    if code != 0:
        params.useRandomCoords = True
        params.randomSeed = int(seed + 1000003)
        code = int(AllChem.EmbedMolecule(molecule, params))
    if code != 0:
        raise RuntimeError("ETKDG embedding failed")
    properties = AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant="MMFF94s")
    if properties is not None:
        force = AllChem.MMFFGetMoleculeForceField(molecule, properties, confId=0)
        force_name = "MMFF94s"
    else:
        force = AllChem.UFFGetMoleculeForceField(molecule, confId=0)
        force_name = "UFF"
    if force is None:
        raise RuntimeError("no force field")
    force.Initialize()
    status = int(force.Minimize(maxIts=500, forceTol=1e-4, energyTol=1e-6))
    energy = float(force.CalcEnergy())
    output = Chem.RemoveHs(molecule)
    output.SetProp("_Name", "QUERY")
    return output, {"embed_code": code, "force_field": force_name,
                    "minimize_status": status, "energy": energy}


def convert_sdf_to_mol2(source: pathlib.Path, destination: pathlib.Path) -> dict:
    command = [str(OBABEL), "-isdf", str(source), "-omol2", "-O", str(destination), "-f", "1", "-l", "1", "-d"]
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=20)
    if completed.returncode != 0 or not destination.is_file() or destination.stat().st_size == 0:
        raise RuntimeError({"obabel_returncode": completed.returncode,
                            "stdout": completed.stdout[-1000:], "stderr": completed.stderr[-1000:]})
    return {"command": command, "stdout": completed.stdout, "stderr": completed.stderr,
            "mol2_sha256": sha256(destination)}


def parse_metrics(stdout: str) -> dict:
    candidates = []
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) == 10:
            try:
                values = [float(x) for x in parts[2:8]]
                sizes = [int(x) for x in parts[8:10]]
            except ValueError:
                continue
            candidates.append((parts, values, sizes))
    if len(candidates) != 1:
        raise RuntimeError("could not uniquely parse LS-align metric row")
    parts, values, sizes = candidates[0]
    return {"query_name": parts[0], "template_name": parts[1],
            "pc_score_query": values[0], "pc_score_template": values[1],
            "pval_query": values[2], "pval_template": values[3],
            "jaccard": values[4], "ls_rmsd_A": values[5],
            "query_size": sizes[0], "template_size": sizes[1]}


def pdb_query_coordinates(path: pathlib.Path, query: Chem.Mol) -> np.ndarray:
    coords, elements = [], []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("TER"):
                break
            if not line.startswith(("ATOM  ", "HETATM")) or line[17:20].strip() != "QUE":
                continue
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            elements.append(line[76:78].strip().capitalize())
    expected = [atom.GetSymbol() for atom in query.GetAtoms()]
    if elements != expected or len(coords) != query.GetNumAtoms():
        raise RuntimeError({"pdb_elements": elements, "expected_elements": expected})
    xyz = np.asarray(coords, dtype=float)
    if not np.isfinite(xyz).all():
        raise RuntimeError("nonfinite LS-align coordinates")
    return xyz


def run_lsalign(stage: pathlib.Path, query_mol2: pathlib.Path, template_mol2: pathlib.Path,
                query: Chem.Mol, mode_name: str, rf: int) -> tuple[Chem.Mol, dict]:
    sup = stage / f"{mode_name}.pdb"
    matrix = stage / f"{mode_name}.matrix.txt"
    command = [str(LSALIGN), str(query_mol2), str(template_mol2), "-rf", str(rf),
               "-md", "0", "-acc", "1", "-oq", "1", "-ot", "1",
               "-o", str(sup), "-m", str(matrix)]
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=60)
    if completed.returncode != 1 or not sup.is_file() or not matrix.is_file():
        raise RuntimeError({"mode": mode_name, "returncode": completed.returncode,
                            "stdout": completed.stdout[-1500:], "stderr": completed.stderr[-1500:]})
    metrics = parse_metrics(completed.stdout)
    xyz = pdb_query_coordinates(sup, query)
    output = Chem.Mol(query)
    conformer = output.GetConformer()
    for atom_index, point in enumerate(xyz):
        conformer.SetAtomPosition(atom_index, Point3D(*map(float, point)))
    if graph_identity(output) != graph_identity(query):
        raise RuntimeError("query graph changed")
    sdf = stage / f"{mode_name}.sdf"
    Chem.MolToMolFile(output, str(sdf))
    return output, {**metrics, "mode": mode_name, "rf": rf, "command": command,
                    "upstream_returncode": completed.returncode,
                    "pose_sha256": sha256(sdf), "superposition_sha256": sha256(sup),
                    "matrix_sha256": sha256(matrix)}


def expected_for(index: int) -> tuple[dict, dict, int, pathlib.Path, pathlib.Path]:
    table = tasks()
    if not 0 <= index < len(table):
        raise IndexError(index)
    _, item, seed = table[index]
    row = item["row"]
    template = DATASET / row["family"] / row["template_id"] / "ligand.sdf"
    protein = DATASET / row["family"] / row["template_id"] / "protein.pdb"
    for path in (MANIFEST, PREREG, template, protein, LSALIGN):
        scope(path)
    if sha256(LSALIGN) != EXPECTED_LSALIGN_SHA256:
        raise RuntimeError("external executable identity drift")
    spec = {
        "task_index": index, "pair": item["pair"], "seed": seed,
        "pair_index": item["pair_index"], "hard_low_2d": False,
        "hard_flex": bool(item["hard_flex"]), "hard_scaffold": bool(item["hard_scaffold"]),
        "query_input_smiles": row["query_input_smiles"],
        "template_path": str(template), "template_sha256": sha256(template),
        "protein_path_recorded_but_not_read": str(protein),
        "manifest_sha256": EXPECTED_MANIFEST_SHA256, "preregister_sha256": sha256(PREREG),
        "lsalign_sha256": EXPECTED_LSALIGN_SHA256, "obabel_sha256": sha256(OBABEL),
        "rdkit_version": rdBase.rdkitVersion, "candidate_modes": [name for name, _ in MODES],
        "selection_rule": "max min(PC-scoreQ,PC-scoreT), then max Jaccard, then min LS-RMSD, then fixed mode name",
        "selector_reads_native": False, "query_native_ligand_read": False,
        "template_ligand_read": True, "template_receptor_read": False,
        "threads_max": 1, "cpu_affinity": [CPU], "cpu_reserve_min": 30,
        "gpu_hidden": True, "timeout_seconds_per_external_call": 60,
    }
    return spec, row, seed, template, protein


def valid_result(path: pathlib.Path, spec: dict) -> bool:
    if not path.is_file():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return payload.get("status") == "success" and all(payload.get(key) == value for key, value in spec.items())


def run(index: int, audit_only: bool = False, fault_after_candidate: bool = False,
        hold_after_candidate: float = 0.0) -> dict:
    global ACTIVE
    if pathlib.Path.cwd().resolve() != ROOT or len(os.sched_getaffinity(0)) < 1:
        raise RuntimeError("scope or CPU reserve precondition failed")
    spec, row, seed, template_path, _protein = expected_for(index)
    if audit_only:
        return {**spec, "status": "audit_pass"}
    os.sched_setaffinity(0, [CPU])
    resource.setrlimit(resource.RLIMIT_AS, (4 * 1024**3, 4 * 1024**3))
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or os.getpid() != os.getpgrp():
        raise RuntimeError("GPU or process-group isolation missing")
    task_root = OUTPUT / spec["pair"] / f"seed_{seed}"
    result_path = task_root / "metadata.json"
    if valid_result(result_path, spec) and (task_root / "selected.sdf").is_file():
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        if payload["selected_pose_sha256"] != sha256(task_root / "selected.sdf"):
            raise RuntimeError("resumed pose checksum drift")
        return {**payload, "status": "resumed"}
    task_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.", dir=task_root.parent))
    ACTIVE = stage
    try:
        query, embed = embed_query(row["query_input_smiles"], seed)
        template = load_single(template_path)
        template.SetProp("_Name", "TEMPLATE")
        query_sdf = stage / "query_initial.sdf"
        template_sdf = stage / "template.sdf"
        Chem.MolToMolFile(query, str(query_sdf))
        Chem.MolToMolFile(template, str(template_sdf))
        query_mol2 = stage / "query.mol2"
        template_mol2 = stage / "template.mol2"
        query_conversion = convert_sdf_to_mol2(query_sdf, query_mol2)
        template_conversion = convert_sdf_to_mol2(template_sdf, template_mol2)
        candidates = []
        for mode_name, rf in MODES:
            molecule, info = run_lsalign(stage, query_mol2, template_mol2, query, mode_name, rf)
            candidates.append((molecule, info))
        selected_molecule, selected = min(
            candidates,
            key=lambda item: (-min(item[1]["pc_score_query"], item[1]["pc_score_template"]),
                              -item[1]["jaccard"], item[1]["ls_rmsd_A"], item[1]["mode"]),
        )
        selected_path = stage / "selected.sdf"
        Chem.MolToMolFile(selected_molecule, str(selected_path))
        ready = {"task_index": index, "pair": spec["pair"], "seed": seed,
                 "selected_mode": selected["mode"], "selected_pose_sha256": sha256(selected_path),
                 "cpu_affinity": sorted(os.sched_getaffinity(0)), "gpu_hidden": True}
        atomic_json(OUTPUT / f"ready_{index}_{os.getpid()}.json", ready)
        if fault_after_candidate:
            raise SystemExit(97)
        if hold_after_candidate > 0:
            time.sleep(hold_after_candidate)
        payload = {**spec, "status": "success", "embedding": embed,
                   "query_initial_sha256": sha256(query_sdf),
                   "query_conversion": query_conversion, "template_conversion": template_conversion,
                   "candidates": [info for _, info in candidates],
                   "selected_mode": selected["mode"],
                   "selected_pose_sha256": sha256(selected_path),
                   "selected_metrics": selected,
                   "SOTA": False}
        atomic_json(stage / "metadata.json", payload)
        if task_root.exists():
            if valid_result(result_path, spec):
                return {**json.loads(result_path.read_text()), "status": "resumed"}
            raise RuntimeError("existing task identity drift")
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
    parser.add_argument("--fault-after-candidate", action="store_true")
    parser.add_argument("--hold-after-candidate", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(run(args.index, args.audit_only, args.fault_after_candidate,
                         args.hold_after_candidate), sort_keys=True))


if __name__ == "__main__":
    main()
