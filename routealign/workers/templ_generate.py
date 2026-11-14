#!/usr/bin/env python3
"""Same-row AlignDockBench rerun of the unmodified TEMPL author core."""

from __future__ import annotations
from routealign import settings

import argparse
import csv
import hashlib
import importlib
import importlib.metadata
import json
import os
import pathlib
import shutil
import sys
import tempfile
import time
import types

from rdkit import Chem


ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
DATASET = ROOT / "data/AlignDockBench"
BUNDLE = ROOT / "tools/templ_source.json"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ADB_TEMPL_OUTPUT",
    str(ROOT / "runs/templ_generate"),
))
EXPECTED_MANIFEST_SHA256 = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
AUTHOR_COMMIT = "727fd306e4066aa39b05a6e51e7833c86d90df98"
AUTHOR_FILES = {
    "mcs.py": "caf216b39a87ce928cf61921269720ac62103bd7791fc3f2ec4aead1e41ee109",
    "scoring.py": "d538ecc50fc4d471a3388eb63e32ccc547ca1a0ff4e1827a5748af6877ce3d76",
    "chemistry.py": "df9f6bd0c3e295349984fa78bb79f92c5510d776b251ba70096163b08f034dd4",
}


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def manifest_rows() -> list[dict[str, str]]:
    if sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError("manifest identity mismatch")
    with MANIFEST.open(newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r["final_primary_eligible"] == "True"]
    rows.sort(key=lambda r: (r["family"], r["query_id"]))
    if len(rows) != 44:
        raise RuntimeError("frozen cohort mismatch")
    return rows


def author_modules(temp_root: pathlib.Path):
    bundle = json.loads(BUNDLE.read_text(encoding="utf-8"))
    if bundle.get("commit") != AUTHOR_COMMIT:
        raise RuntimeError("TEMPL author commit mismatch")
    package_name = "routealign_templ_author_core"
    package = temp_root / package_name
    package.mkdir()
    (package / "__init__.py").write_text("", encoding="utf-8")
    for name, expected in AUTHOR_FILES.items():
        raw = bundle["contents"]["templ_pipeline/core/" + name]
        if hashlib.sha256(raw.encode()).hexdigest() != expected:
            raise RuntimeError(f"TEMPL author source mismatch: {name}")
        if name == "mcs.py":
            if raw.count("ps.coordMap = coordMap") != 1 or raw.count("ps.coordMap = cmap") != 1:
                raise RuntimeError("unexpected author CoordMap source")
            raw = raw.replace("ps.coordMap = coordMap", "ps.SetCoordMap(coordMap)")
            raw = raw.replace("ps.coordMap = cmap", "ps.SetCoordMap(cmap)")
        (package / name).write_text(raw, encoding="utf-8")
    sys.path.insert(0, str(temp_root))
    try:
        mcs = importlib.import_module(package_name + ".mcs")
        scoring = importlib.import_module(package_name + ".scoring")
    finally:
        sys.path.remove(str(temp_root))
    return mcs, scoring


def canonical(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)


def load_template(path: pathlib.Path) -> Chem.Mol:
    values = [m for m in Chem.SDMolSupplier(str(path), removeHs=True) if m is not None]
    if len(values) != 1 or values[0].GetNumConformers() != 1:
        raise RuntimeError(f"invalid template SDF: {path}")
    return values[0]


def valid_existing(path: pathlib.Path, expected: dict) -> bool:
    if not path.exists():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") != "success" or any(payload.get(k) != v for k, v in expected.items()):
            return False
        return all(sha256(path.parent / name) == value for name, value in payload["output_sha256"].items())
    except Exception:
        return False


def run_task(index: int, fault_after_generation: bool, hold_after_generation: float) -> dict:
    if pathlib.Path.cwd().resolve() != ROOT or (os.cpu_count() or 0) < 1:
        raise RuntimeError("repository or CPU reserve mismatch")
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise RuntimeError("TEMPL baseline must hide GPUs")
    rows = manifest_rows()
    if not 0 <= index < len(rows):
        raise IndexError(index)
    row = rows[index]
    task_root = OUTPUT / f"{row['family']}__{row['query_id']}"
    result_path = task_root / "result.json"
    template_path = DATASET / row["family"] / row["template_id"] / "ligand.sdf"
    expected = {
        "task_index": index,
        "family": row["family"],
        "query_id": row["query_id"],
        "manifest_sha256": EXPECTED_MANIFEST_SHA256,
        "template_sha256": sha256(template_path),
        "worker_sha256": sha256(pathlib.Path(__file__)),
        "author_commit": AUTHOR_COMMIT,
        "rdkit_version": importlib.metadata.version("rdkit"),
        "numpy_version": importlib.metadata.version("numpy"),
    }
    if valid_existing(result_path, expected):
        return {"status": "resumed", "output": str(task_root), **expected}
    if task_root.exists():
        raise RuntimeError(f"incomplete or drifting output exists: {task_root}")

    template = load_template(template_path)
    query = Chem.MolFromSmiles(row["query_input_smiles"])
    if query is None or canonical(query) == canonical(template):
        raise RuntimeError("query graph invalid or identical to template")
    temp_module_root = pathlib.Path(tempfile.mkdtemp(prefix="routealign_templ_module_", dir="/tmp"))
    try:
        mcs, scoring = author_modules(temp_module_root)
        target = Chem.AddHs(query)
        _, smarts, details = mcs.find_mcs(target, [template], return_details=True)
        generated = mcs.constrained_embed(
            target, template, smarts, n_conformers=50,
            n_workers_pipeline=1, enable_optimization=False,
        )
        if generated is None or generated.GetNumConformers() < 1:
            raise RuntimeError("TEMPL author constrained embedding failed")
        ranked = scoring.select_best(
            generated, template, n_workers=1,
            return_all_ranked=True, align_metric="combo",
        )
        if not ranked:
            raise RuntimeError("TEMPL author ranking failed")
        blocks = []
        scores = []
        for rank, (conformer_id, values, mol) in enumerate(ranked):
            output = Chem.RemoveHs(Chem.Mol(mol))
            if canonical(output) != canonical(query) or output.GetNumConformers() != 1:
                raise RuntimeError("TEMPL changed the complete query graph")
            output.SetProp("_Name", f"{row['family']}__{row['query_id']}_TEMPL_rank_{rank:03d}")
            blocks.append((Chem.MolToMolBlock(output) + "\n$$$$\n").encode())
            scores.append({"rank": rank, "source_conformer_id": int(conformer_id), "scores": {k: float(v) for k, v in values.items()}})
    finally:
        for key in [k for k in list(sys.modules) if k == "routealign_templ_author_core" or k.startswith("routealign_templ_author_core.")]:
            del sys.modules[key]
        shutil.rmtree(temp_module_root, ignore_errors=True)

    if fault_after_generation:
        raise SystemExit(97)
    if hold_after_generation:
        time.sleep(hold_after_generation)

    task_root.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=task_root.name + ".staging.", dir=task_root.parent))
    try:
        (stage / "selected.sdf").write_bytes(blocks[0])
        (stage / "ranked.sdf").write_bytes(b"".join(blocks))
        payload = {
            **expected,
            "status": "success",
            "method": "TEMPL author core, single supplied template, combo top-1",
            "author_core_sha256": AUTHOR_FILES,
            "coordmap_api_fix": "two removed RDKit attributes replaced by documented SetCoordMap calls",
            "author_seed": 42,
            "generated_conformers": int(generated.GetNumConformers()),
            "ranked_conformers": len(ranked),
            "mcs_details": details,
            "ranked_scores": scores,
            "n_conformers": 50,
            "enable_optimization": False,
            "query_native_ligand_read": False,
            "target_pose_metric_read": False,
            "rdkit_version": importlib.metadata.version("rdkit"),
            "numpy_version": importlib.metadata.version("numpy"),
            "output_sha256": {
                "selected.sdf": sha256(stage / "selected.sdf"),
                "ranked.sdf": sha256(stage / "ranked.sdf"),
            },
        }
        (stage / "result.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        stage.replace(task_root)
    finally:
        if stage.exists():
            shutil.rmtree(stage, ignore_errors=True)
    return {"status": "computed", "output": str(task_root), **expected}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--list", action="store_true")
    p.add_argument("--index", type=int)
    p.add_argument("--fault-after-generation", action="store_true")
    p.add_argument("--hold-after-generation", type=float, default=0.0)
    a = p.parse_args()
    if a.list:
        for i, row in enumerate(manifest_rows()):
            print(i, row["family"], row["query_id"])
        return
    if a.index is None:
        p.error("--index is required")
    print(json.dumps(run_task(a.index, a.fault_after_generation, a.hold_after_generation), sort_keys=True))


if __name__ == "__main__":
    main()
