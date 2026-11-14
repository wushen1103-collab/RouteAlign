#!/usr/bin/env python3
"""Outcome-blind hard-core grafting under the AlignDockBench cross-docking protocol."""

from __future__ import annotations
from routealign import settings

import argparse
import csv
import hashlib
import json
import math
import os
import pathlib
import tempfile
import time

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdFMCS, rdShapeHelpers
from rdkit.Geometry import Point3D


ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
DATASET = ROOT / "data/AlignDockBench"
OUTPUT = pathlib.Path(os.environ.get("ROUTEALIGN_ADB_GRAFT_OUTPUT", str(ROOT / "runs/graft")))
EXPECTED_MANIFEST_SHA256 = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
SEEDS = (101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163,
         167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229)
VARIANTS = ("energy", "shape", "hybrid", "pocket_hybrid")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: pathlib.Path, payload: dict) -> None:
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


def manifest_rows() -> list[dict[str, str]]:
    if sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError("manifest identity mismatch")
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["final_primary_eligible"] == "True"]
    rows.sort(key=lambda row: (row["family"], row["query_id"]))
    if len(rows) != 44:
        raise RuntimeError("cohort mismatch")
    return rows


def tasks() -> list[tuple[int, dict[str, str], int]]:
    return [(i, row, seed) for i, (row, seed) in enumerate((r, s) for r in manifest_rows() for s in SEEDS)]


def load_single(path: pathlib.Path) -> Chem.Mol:
    molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if mol is not None]
    if len(molecules) != 1 or molecules[0].GetNumConformers() != 1:
        raise RuntimeError(f"invalid SDF: {path}")
    return molecules[0]


def protein_heavy_xyz(path: pathlib.Path) -> np.ndarray:
    points = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if element == "H":
                continue
            try:
                points.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    if not points:
        raise RuntimeError(f"no protein heavy atoms: {path}")
    return np.asarray(points, dtype=float)


def optimize_fixed_core(mol: Chem.Mol, fixed: tuple[int, ...]) -> tuple[float, str, int]:
    properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if properties is not None:
        force_field = AllChem.MMFFGetMoleculeForceField(mol, properties, confId=0)
        force_name = "MMFF94s"
    else:
        force_field = AllChem.UFFGetMoleculeForceField(mol, confId=0)
        force_name = "UFF"
    if force_field is None:
        raise RuntimeError("no force field")
    for atom_index in fixed:
        force_field.AddFixedPoint(int(atom_index))
    force_field.Initialize()
    status = int(force_field.Minimize(maxIts=500, forceTol=1e-4, energyTol=1e-6))
    return float(force_field.CalcEnergy()), force_name, status


def clash_penalty(mol: Chem.Mol, protein_xyz: np.ndarray) -> float:
    xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    minimum = np.full(len(xyz), np.inf)
    for start in range(0, len(protein_xyz), 2048):
        block = protein_xyz[start : start + 2048]
        distances = np.linalg.norm(xyz[:, None, :] - block[None, :, :], axis=2)
        minimum = np.minimum(minimum, distances.min(axis=1))
    return float(np.mean(np.maximum(0.0, 1.8 - minimum) ** 2))


def graph_identity(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def build_candidates(row: dict[str, str], seed: int) -> tuple[list[dict], Chem.Mol]:
    template_path = DATASET / row["family"] / row["template_id"] / "ligand.sdf"
    # Cross-docking may use the protein co-crystallized with the supplied template only.
    protein_path = DATASET / row["family"] / row["template_id"] / "protein.pdb"
    template = load_single(template_path)
    query = Chem.MolFromSmiles(row["query_input_smiles"])
    if query is None:
        raise RuntimeError("query graph parse failed")
    mcs = rdFMCS.FindMCS(
        [query, template], timeout=60, atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact, ringMatchesRingOnly=True,
        completeRingsOnly=True, matchValences=True,
    )
    if mcs.canceled or mcs.numAtoms != int(row["mcs_atoms"]):
        raise RuntimeError("frozen MCS mismatch")
    pattern = Chem.MolFromSmarts(mcs.smartsString)
    query_matches = sorted(query.GetSubstructMatches(pattern, uniquify=False, maxMatches=128, useChirality=True))
    template_matches = sorted(template.GetSubstructMatches(pattern, uniquify=False, maxMatches=128, useChirality=True))
    mapping_pairs = [(q, t) for q in query_matches for t in template_matches]
    if not mapping_pairs:
        raise RuntimeError("no stereo-compatible MCS mapping")
    mapping_pairs = mapping_pairs[:32]
    template_xyz = np.asarray(template.GetConformer().GetPositions(), dtype=float)
    pocket_xyz = protein_heavy_xyz(protein_path)
    candidates = []
    diagnostics = {"embed_failed": 0, "force_field_failed": 0, "graph_changed": 0, "nonfinite": 0}
    force_field_errors = []
    for mapping_index, (query_match, template_match) in enumerate(mapping_pairs):
        attempt_schedule = (
            (int(seed + 1009 * mapping_index), True),
            (int(seed + 1009 * mapping_index + 1000), True),
            (int(seed + 1009 * mapping_index + 2000), True),
            (int(seed + 1009 * mapping_index + 3000), False),
        )
        for conformer_index, (attempt_seed, use_random_coords) in enumerate(attempt_schedule):
            molecule = Chem.AddHs(query)
            parameters = AllChem.ETKDGv3()
            parameters.randomSeed = attempt_seed
            parameters.useRandomCoords = use_random_coords
            parameters.numThreads = 1
            if hasattr(parameters, "maxIterations"):
                parameters.maxIterations = 1000
            parameters.SetCoordMap(
                {int(q): Point3D(*map(float, template_xyz[t])) for q, t in zip(query_match, template_match)}
            )
            code = int(AllChem.EmbedMolecule(molecule, parameters))
            if code != 0:
                diagnostics["embed_failed"] += 1
                continue
            conformer = molecule.GetConformer()
            for query_atom, template_atom in zip(query_match, template_match):
                conformer.SetAtomPosition(int(query_atom), Point3D(*map(float, template_xyz[template_atom])))
            try:
                energy, force_field, minimize_status = optimize_fixed_core(molecule, tuple(query_match))
            except (RuntimeError, ValueError) as exc:
                diagnostics["force_field_failed"] += 1
                if len(force_field_errors) < 5:
                    force_field_errors.append(repr(exc))
                continue
            output = Chem.RemoveHs(molecule)
            if graph_identity(output) != graph_identity(query):
                diagnostics["graph_changed"] += 1
                continue
            output_xyz = np.asarray(output.GetConformer().GetPositions(), dtype=float)
            if not np.isfinite(output_xyz).all():
                diagnostics["nonfinite"] += 1
                continue
            core_deviation = float(
                np.sqrt(np.mean(np.sum((output_xyz[list(query_match)] - template_xyz[list(template_match)]) ** 2, axis=1)))
            )
            candidate = {
                "mol": output,
                "mapping_index": mapping_index,
                "conformer_index": conformer_index,
                "attempt_seed": attempt_seed,
                "use_random_coords": use_random_coords,
                "query_match": list(query_match),
                "template_match": list(template_match),
                "force_field": force_field,
                "minimize_status": minimize_status,
                "energy": energy,
                "energy_per_heavy_atom": energy / output.GetNumHeavyAtoms(),
                "shape_distance": float(rdShapeHelpers.ShapeTanimotoDist(output, template)),
                "protein_clash_penalty": clash_penalty(output, pocket_xyz),
                "core_coordinate_deviation_A": core_deviation,
            }
            candidates.append(candidate)
    if not candidates:
        raise RuntimeError(
            "all candidate constructions failed: "
            + json.dumps({"mapping_pairs": len(mapping_pairs), "diagnostics": diagnostics, "force_field_errors": force_field_errors})
        )
    return candidates, template


def rank_candidates(candidates: list[dict]) -> dict[str, int]:
    energies = np.asarray([candidate["energy_per_heavy_atom"] for candidate in candidates])
    energy_delta = energies - energies.min()
    for candidate, delta in zip(candidates, energy_delta):
        candidate["hybrid_score"] = float(candidate["shape_distance"] + 0.02 * delta)
        candidate["pocket_hybrid_score"] = float(
            candidate["shape_distance"] + 0.02 * delta + 0.5 * candidate["protein_clash_penalty"]
        )
    criteria = {
        "energy": "energy_per_heavy_atom",
        "shape": "shape_distance",
        "hybrid": "hybrid_score",
        "pocket_hybrid": "pocket_hybrid_score",
    }
    return {
        variant: min(range(len(candidates)), key=lambda i: (candidates[i][field], candidates[i]["mapping_index"], candidates[i]["conformer_index"]))
        for variant, field in criteria.items()
    }


def valid_existing(path: pathlib.Path, expected: dict) -> bool:
    if not path.is_file():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return payload.get("status") == "success" and all(payload.get(key) == value for key, value in expected.items())


def run_task(index: int, fault_after_candidates: bool, hold_after_candidates: float) -> dict:
    table = tasks()
    if not 0 <= index < len(table):
        raise IndexError(index)
    _, row, seed = table[index]
    task_root = OUTPUT / f"{row['family']}__{row['query_id']}" / f"seed_{seed}"
    metadata_path = task_root / "metadata.json"
    template_path = DATASET / row["family"] / row["template_id"] / "ligand.sdf"
    protein_path = DATASET / row["family"] / row["template_id"] / "protein.pdb"
    expected = {
        "task_index": index, "family": row["family"], "query_id": row["query_id"], "seed": seed,
        "manifest_sha256": EXPECTED_MANIFEST_SHA256,
        "template_sha256": sha256(template_path), "protein_sha256": sha256(protein_path),
    }
    if valid_existing(metadata_path, expected):
        return {"status": "resumed", "output": str(task_root), **expected}
    if task_root.exists():
        raise RuntimeError(f"incomplete output exists and will not be overwritten: {task_root}")
    candidates, _ = build_candidates(row, seed)
    selection = rank_candidates(candidates)
    if fault_after_candidates:
        raise SystemExit(97)
    if hold_after_candidates > 0:
        time.sleep(hold_after_candidates)
    task_root.parent.mkdir(parents=True, exist_ok=True)
    staging = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.", dir=task_root.parent))
    try:
        selected = {}
        for variant, candidate_index in selection.items():
            path = staging / f"pose_{variant}.sdf"
            writer = Chem.SDWriter(str(path))
            writer.write(candidates[candidate_index]["mol"])
            writer.close()
            selected[variant] = {
                "candidate_index": candidate_index,
                "pose_sha256": sha256(path),
                **{key: value for key, value in candidates[candidate_index].items() if key != "mol"},
            }
        payload = {
            **expected, "status": "success", "method": "hard MCS coordinate graft + fixed-core force-field relaxation + template-pocket reranking",
            "variants": list(VARIANTS), "candidate_count": len(candidates), "selected": selected,
            "protein_input_role": "template co-crystal protein (cross-docking)",
            "query_native_protein_read": False, "query_native_ligand_read": False,
            "target_pose_metric_read": False, "rdkit_version": rdBase.rdkitVersion,
        }
        atomic_json(staging / "metadata.json", payload)
        os.replace(staging, task_root)
    except BaseException:
        raise
    return {"status": "computed", "output": str(task_root), **expected}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--fault-after-candidates", action="store_true")
    parser.add_argument("--hold-after-candidates", type=float, default=0.0)
    args = parser.parse_args()
    if args.list:
        for index, row, seed in tasks():
            print(index, row["family"], row["query_id"], seed)
        return
    if args.index is None:
        parser.error("--index is required")
    print(json.dumps(run_task(args.index, args.fault_after_candidates, args.hold_after_candidates), sort_keys=True))


if __name__ == "__main__":
    main()

