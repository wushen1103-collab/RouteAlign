#!/usr/bin/env python3
"""Full and confirmatory same-row analysis for the frozen 44-pair method."""
from __future__ import annotations
from routealign import settings

import argparse
import csv
import hashlib
import itertools
import json
import os
import pathlib
import shutil
import signal
import tempfile
import time

import numpy as np


ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
PREREG = ROOT / "config/protocol.json"
ENDPOINT = ROOT / "runs/medoid_evaluate"
TEMPL = ROOT / "runs/templ_evaluate"
IKTOS = ROOT / "data/published_baselines"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ADB_GROUP_ANALYSIS_OUTPUT",
    str(ROOT / "runs/uniform_statistics"),
)).resolve()
EXPECTED_MANIFEST = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
EXPECTED_ALIGN = "1062dbf414525b0d0093ce59c4a85ddeefbf338ad8ead8c7c0c1e12d99064be6"
EXPECTED_DOCK = "fcd080d8cfebe1801919f0ec5a4fda188cc08dfeb3bfeeaa0d29136e2b2ec0b3"
BOOTSTRAP_SEED = 20261027
BOOTSTRAP_REPLICATES = 20000


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def csv_rows(path):
    with pathlib.Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sample_std(values):
    return float(np.std(np.asarray(values, dtype=float), ddof=1))


def exact_family_signflip(rows, pair_delta):
    families = sorted({row["family"] for row in rows})
    values = np.asarray([
        np.mean([pair_delta[i] for i, row in enumerate(rows) if row["family"] == family])
        for family in families
    ], dtype=float)
    observed = float(values.mean())
    if len(values) <= 16:
        count = 0
        total = 0
        for signs in itertools.product((-1.0, 1.0), repeat=len(values)):
            statistic = float(np.mean(values * np.asarray(signs)))
            count += statistic <= observed + 1e-15
            total += 1
        return {"family_count": len(families), "observed_family_mean_delta_A": observed,
                "alternative": "proposed RMSD lower", "one_sided_p": (count + 1) / (total + 1),
                "method": "exact", "sign_patterns": total}
    draws = 200000
    rng = np.random.default_rng(BOOTSTRAP_SEED + len(values))
    signs = rng.choice(np.asarray((-1.0, 1.0)), size=(draws, len(values)), replace=True)
    statistics = signs @ values / len(values)
    count = int(np.sum(statistics <= observed + 1e-15))
    return {"family_count": len(families), "observed_family_mean_delta_A": observed,
            "alternative": "proposed RMSD lower", "one_sided_p": (count + 1) / (draws + 1),
            "method": "fixed-seed Monte Carlo", "draws": draws,
            "random_seed": BOOTSTRAP_SEED + len(values)}


def bootstrap(rows, proposed, baseline):
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    families = sorted({row["family"] for row in rows})
    family_indices = {family: np.asarray([i for i, row in enumerate(rows) if row["family"] == family])
                      for family in families}
    dr = np.empty(BOOTSTRAP_REPLICATES)
    ds = np.empty(BOOTSTRAP_REPLICATES)
    for iteration in range(BOOTSTRAP_REPLICATES):
        sampled_families = rng.choice(families, len(families), replace=True)
        indices = np.concatenate([family_indices[family] for family in sampled_families])
        groups = rng.integers(0, proposed.shape[1], proposed.shape[1])
        p = proposed[indices][:, groups]
        b = baseline[indices]
        dr[iteration] = p.mean() - b.mean()
        ds[iteration] = (p <= 2.0).mean() - (b <= 2.0).mean()
    return {"replicates": BOOTSTRAP_REPLICATES, "random_seed": BOOTSTRAP_SEED,
            "cluster_unit": "protein family", "resamples_independent_seed_groups": True,
            "delta_mean_rmsd_A_ci95": [float(np.quantile(dr, .025)), float(np.quantile(dr, .975))],
            "probability_delta_mean_rmsd_ge_0": float(np.mean(dr >= 0)),
            "delta_success_ci95": [float(np.quantile(ds, .025)), float(np.quantile(ds, .975))]}


def load():
    identities = {
        "manifest": (MANIFEST, EXPECTED_MANIFEST),
        "iktos_align": (IKTOS / "3DAlign_benchmark_results.csv", EXPECTED_ALIGN),
        "iktos_dock": (IKTOS / "3DDock_benchmark_results.csv", EXPECTED_DOCK),
    }
    for label, (path, expected) in identities.items():
        if sha256(path) != expected:
            raise RuntimeError(label + " identity drift")
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    if not prereg.get("stagewise_frozen_before_full_expansion_batch") or len(prereg["pairs"]) != 44:
        raise RuntimeError("preregistration mismatch")
    manifest = [row for row in csv_rows(MANIFEST) if row["final_primary_eligible"] == "True"]
    by_pair = {f"{row['family']}__{row['query_id']}": row for row in manifest}
    align = {(row["batch_job_id"].lower(), row["ligand_id"].lower()): row["rmsd"]
             for row in csv_rows(identities["iktos_align"][0])}
    dock = {(row["batch_job_id"].lower(), row["ligand_id"].lower()): row["rmsd"]
            for row in csv_rows(identities["iktos_dock"][0])}
    records = []
    input_hashes = [EXPECTED_MANIFEST, EXPECTED_ALIGN, EXPECTED_DOCK, sha256(PREREG)]
    for pair_spec in prereg["pairs"]:
        pair = pair_spec["pair"]
        row = by_pair[pair]
        values = []
        pb = []
        joint = []
        for group in range(1, 6):
            path = ENDPOINT / pair / f"group_{group}" / "result.json"
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("status") != "success" or payload.get("pair") != pair or payload.get("group_index") != group:
                raise RuntimeError("endpoint identity mismatch: " + str(path))
            values.append(float(payload["rmsd_A"]))
            pb.append(bool(payload["PB20_pass"]))
            joint.append(bool(payload["joint_success"]))
            input_hashes.append(sha256(path))
        key = (row["family"].lower(), row["query_id"].lower())
        templ_path = TEMPL / pair / "result.json"
        templ = json.loads(templ_path.read_text(encoding="utf-8"))
        input_hashes.append(sha256(templ_path))
        records.append({**row, **pair_spec, "groups": values, "pb": pb, "joint": joint,
                        "baselines": {"Iktos-3DAlign": float(align[key]), "Iktos-3DDock": float(dock[key]),
                                      "TEMPL-author": float(templ["rmsd_A"])}})
    digest = hashlib.sha256("".join(sorted(input_hashes)).encode("ascii")).hexdigest()
    return records, digest


def holm_adjust(comparisons):
    ordered = sorted(comparisons, key=lambda item: item[2])
    running = 0.0
    adjusted = {}
    total = len(ordered)
    for rank, (stratum, method, pvalue) in enumerate(ordered):
        running = max(running, min(1.0, pvalue * (total - rank)))
        adjusted[(stratum, method)] = running
    return adjusted


def analyze():
    records, digest = load()
    predicates = {
        "primary_all44": lambda row: True,
        "confirmatory_non_low_2d": lambda row: not bool(row["hard_low_2d"]),
        "hard_low_2d_similarity": lambda row: bool(row["hard_low_2d"]),
        "hard_scaffold_change": lambda row: bool(row["hard_scaffold"]),
        "hard_upper_flexibility_band": lambda row: bool(row["hard_flex"]),
        "hard_low_2d_and_flexible": lambda row: bool(row["hard_low_2d"] and row["hard_flex"]),
        "hard_low_2d_scaffold_and_flexible": lambda row: bool(row["hard_low_2d"] and row["hard_flex"] and row["hard_scaffold"]),
    }
    result = {"analysis_version": "group-medoid-full-confirmatory-v5", "input_digest": digest,
              "protocol": {"method_replicates": 5, "replicate_unit": "disjoint five-seed groups",
                           "success_threshold_A": 2.0, "difficulty_strata_overlap": True,
                           "prediction_selector_reads_native": False,
                           "baseline_runs_per_pair": 1,
                           "timeout_fallback": "2ica__3m6f used a preregistered deterministic non-random embedding fallback after the full schedule exceeded 900 seconds",
                           "baseline_sources": {"Iktos-3DAlign": "official same-row public output, commit 747d842",
                                                "Iktos-3DDock": "official same-row public output, commit 747d842",
                                                "TEMPL-author": "self-rerun author code, commit 727fd30, same rows"}},
              "strata": {}}
    raw_p = []
    for name, predicate in predicates.items():
        rows = [row for row in records if predicate(row)]
        proposed = np.asarray([row["groups"] for row in rows], dtype=float)
        pb = np.asarray([row["pb"] for row in rows], dtype=bool)
        joint = np.asarray([row["joint"] for row in rows], dtype=bool)
        group_rmsd = proposed.mean(axis=0)
        group_success = (proposed <= 2.0).mean(axis=0)
        group_pb = pb.mean(axis=0)
        group_joint = joint.mean(axis=0)
        entry = {"n_pairs": len(rows), "n_families": len({row["family"] for row in rows}),
                 "proposed": {"success_mean": float(group_success.mean()), "success_std": sample_std(group_success),
                              "rmsd_A_mean": float(group_rmsd.mean()), "rmsd_A_std": sample_std(group_rmsd),
                              "PB20_mean": float(group_pb.mean()), "PB20_std": sample_std(group_pb),
                              "joint_mean": float(group_joint.mean()), "joint_std": sample_std(group_joint),
                              "per_group_mean_rmsd_A": group_rmsd.tolist()}, "baselines": {}}
        dominance = []
        for method in ("Iktos-3DAlign", "Iktos-3DDock", "TEMPL-author"):
            baseline = np.asarray([row["baselines"][method] for row in rows], dtype=float)
            pair_delta = proposed.mean(axis=1) - baseline
            signflip = exact_family_signflip(rows, pair_delta)
            raw_p.append((name, method, signflip["one_sided_p"]))
            all_groups = [bool(group_success[i] >= (baseline <= 2.0).mean() and group_rmsd[i] < baseline.mean())
                          for i in range(5)]
            entry["baselines"][method] = {
                "success": float((baseline <= 2.0).mean()), "mean_rmsd_A": float(baseline.mean()),
                "delta_proposed_mean_rmsd_A": float(proposed.mean() - baseline.mean()),
                "pair_win_rate": float((proposed.mean(axis=1) < baseline).mean()),
                "all_five_groups_dominate": all(all_groups), "per_group_dominance": all_groups,
                "hierarchical_bootstrap": bootstrap(rows, proposed, baseline),
                "family_signflip": signflip,
            }
            dominance.append(all(all_groups))
        entry["same_row_all_five_group_dominance"] = all(dominance)
        entry["claim_scope"] = "repeated conditional same-row leader" if all(dominance) else "no dominance claim"
        result["strata"][name] = entry
    adjusted = holm_adjust(raw_p)
    for name, method, _ in raw_p:
        result["strata"][name]["baselines"][method]["family_signflip"]["holm_adjusted_p_across_21"] = adjusted[(name, method)]
    return result, digest


def run(args):
    result, digest = analyze()
    if args.audit_only:
        print(json.dumps({"status": "audit_pass", "input_digest": digest,
                          "strata": {key: value["n_pairs"] for key, value in result["strata"].items()},
                          "cpu_threads_max": 1, "gpu_hidden": True}, sort_keys=True))
        return 0
    receipt = OUTPUT / "completion.json"
    worker_hash = sha256(__file__)
    if receipt.is_file():
        old = json.loads(receipt.read_text(encoding="utf-8"))
        if old.get("input_digest") == digest and old.get("worker_sha256") == worker_hash:
            print(json.dumps({"status": "resumed", "output": str(OUTPUT)}, sort_keys=True))
            return 0
        raise RuntimeError("resume identity drift")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    stage = pathlib.Path(tempfile.mkdtemp(prefix=OUTPUT.name + ".staging.", dir=OUTPUT.parent))
    def controlled_stop(_sig, _frame):
        raise TimeoutError("controlled timeout before atomic commit")
    signal.signal(signal.SIGTERM, controlled_stop)
    try:
        (stage / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if args.fault_after_write:
            raise SystemExit(97)
        delay = float(os.environ.get("ROUTEALIGN_ANALYSIS_DELAY_BEFORE_COMMIT", "0"))
        if delay:
            time.sleep(delay)
        (stage / "completion.json").write_text(json.dumps({"status": "complete", "input_digest": digest,
                                                            "worker_sha256": worker_hash}, indent=2, sort_keys=True) + "\n")
        os.replace(stage, OUTPUT)
        print(json.dumps({"status": "complete", "output": str(OUTPUT)}, sort_keys=True))
        return 0
    finally:
        if stage.exists():
            shutil.rmtree(stage)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--fault-after-write", action="store_true")
    raise SystemExit(run(parser.parse_args()))
