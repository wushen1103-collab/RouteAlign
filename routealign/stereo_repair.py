"""Graph-guided correction of coordinate stereochemistry.

The expected molecular graph supplies the stereochemical labels.  Candidate
coordinates are the only geometry read by this module; a native pose is neither
accepted nor required.  The current implementation handles tetrahedral centres
with one implicit hydrogen by reflecting the smallest acyclic substituent across
the plane through the centre and the two retained heavy-atom neighbours.
"""
from __future__ import annotations

import itertools
import numpy as np
from rdkit import Chem


def infer_from_coordinates(mol: Chem.Mol) -> Chem.Mol:
    """Return a copy whose tetrahedral labels are inferred from conformer 0."""
    inferred = Chem.Mol(mol)
    Chem.RemoveStereochemistry(inferred)
    Chem.AssignAtomChiralTagsFromStructure(
        inferred, confId=0, replaceExistingTags=True
    )
    Chem.AssignStereochemistry(inferred, cleanIt=True, force=True)
    return inferred


def isomeric_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def _smallest_acyclic_substituent(mol: Chem.Mol, centre: int):
    neighbours = [
        atom.GetIdx()
        for atom in mol.GetAtomWithIdx(centre).GetNeighbors()
        if atom.GetAtomicNum() > 1
    ]
    if len(neighbours) != 3:
        raise ValueError(
            f"atom {centre} has {len(neighbours)} heavy neighbours; expected 3"
        )
    components = []
    for moving in neighbours:
        seen, stack, component = {centre}, [moving], set()
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            component.add(current)
            stack.extend(
                atom.GetIdx()
                for atom in mol.GetAtomWithIdx(current).GetNeighbors()
                if atom.GetIdx() not in seen
            )
        if not any(index in component for index in neighbours if index != moving):
            components.append((len(component), moving, sorted(component)))
    if not components:
        raise ValueError(f"no acyclic substituent at atom {centre}")
    _, moving, component = min(components)
    return moving, component, [index for index in neighbours if index != moving]


def reflect_substituent(mol: Chem.Mol, centre: int) -> list[int]:
    """Invert one centre while preserving its bond lengths and bond angles."""
    _, component, retained = _smallest_acyclic_substituent(mol, centre)
    conformer = mol.GetConformer()
    origin = np.asarray(conformer.GetAtomPosition(centre), dtype=float)
    points = np.asarray(
        [conformer.GetAtomPosition(index) for index in retained], dtype=float
    )
    normal = np.cross(points[0] - origin, points[1] - origin)
    norm = float(np.linalg.norm(normal))
    if norm <= 1e-8:
        raise ValueError(f"degenerate neighbour plane at atom {centre}")
    normal /= norm
    for index in component:
        point = np.asarray(conformer.GetAtomPosition(index), dtype=float)
        reflected = point - 2.0 * np.dot(point - origin, normal) * normal
        conformer.SetAtomPosition(index, reflected.tolist())
    return component


def repair_candidate(candidate: Chem.Mol, expected: Chem.Mol):
    """Correct candidate coordinates to exactly match ``expected`` stereo.

    Returns the repaired molecule and a compact audit record.  All compatible
    graph maps are considered, which makes the centre discovery independent of
    candidate atom order.  Symmetry-aware candidate-to-candidate RMSD remains a
    separate operation in :mod:`routealign.atom_correspondence`.
    """
    matches = candidate.GetSubstructMatches(
        expected, uniquify=False, useChirality=False, maxMatches=100000
    )
    if not matches:
        raise ValueError("expected graph does not map to candidate")
    expected_centres = [
        index
        for index, _ in Chem.FindMolChiralCenters(
            expected, includeUnassigned=False, includeCIP=True
        )
    ]
    candidate_centres = set()
    for match in matches:
        for expected_index in expected_centres:
            index = match[expected_index]
            atom = candidate.GetAtomWithIdx(index)
            if sum(n.GetAtomicNum() > 1 for n in atom.GetNeighbors()) == 3:
                candidate_centres.add(index)
    expected_label = isomeric_smiles(expected)
    original_xyz = np.asarray(candidate.GetConformer().GetPositions(), dtype=float)
    best = None
    ordered_centres = sorted(candidate_centres)
    for bits in itertools.product([0, 1], repeat=len(ordered_centres)):
        trial = Chem.Mol(candidate)
        flipped = []
        for flag, centre in zip(bits, ordered_centres):
            if flag:
                reflect_substituent(trial, centre)
                flipped.append(centre)
        inferred = infer_from_coordinates(trial)
        if isomeric_smiles(inferred) != expected_label:
            continue
        xyz = np.asarray(trial.GetConformer().GetPositions(), dtype=float)
        displacement = np.sqrt(np.sum((xyz - original_xyz) ** 2, axis=1))
        record = (
            float(np.sqrt(np.mean(displacement**2))),
            float(displacement.max()),
            flipped,
            inferred,
        )
        if best is None or record[:2] < best[:2]:
            best = record
    if best is None:
        raise RuntimeError("no exact local stereochemical correction found")
    movement, maximum, flipped, repaired = best
    return repaired, {
        "expected_isomeric_smiles": expected_label,
        "repaired_isomeric_smiles": isomeric_smiles(repaired),
        "movement_rmsd_A": movement,
        "maximum_atom_movement_A": maximum,
        "flipped_atom_indices": flipped,
        "native_pose_read": False,
    }
