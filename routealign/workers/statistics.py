#!/usr/bin/env python3
"""Audited conditional analysis including ROSHAMBO JCIM 2024 endpoints."""
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
LOW2D_ENDPOINT = ROOT / "runs/medoid_evaluate"
NONLOW2D_ENDPOINT = ROOT / "runs/receptor_evaluate"
TEMPL = ROOT / "runs/templ_evaluate"
IKTOS = ROOT / "data/published_baselines"
LSALIGN_GENERATION = ROOT / "runs/lsalign_generate"
LSALIGN_ENDPOINT = ROOT / "runs/lsalign_evaluate"
LSALIGN_PB = ROOT / "runs/lsalign_validate"
ROSHAMBO_ENDPOINT = ROOT / "runs/roshambo_evaluate"
ROSHAMBO_PB = ROOT / "runs/roshambo_validate"
OUTPUT = pathlib.Path(os.environ.get(
    "ROUTEALIGN_ADB_GROUP_ANALYSIS_OUTPUT",
    str(ROOT / "runs/statistics"),
)).resolve()
EXPECTED_MANIFEST = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
EXPECTED_ALIGN = "1062dbf414525b0d0093ce59c4a85ddeefbf338ad8ead8c7c0c1e12d99064be6"
EXPECTED_DOCK = "fcd080d8cfebe1801919f0ec5a4fda188cc08dfeb3bfeeaa0d29136e2b2ec0b3"
BOOTSTRAP_SEED = 20261027
BOOTSTRAP_REPLICATES = 20000
SEEDS = (101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
         151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
         199, 211, 223, 227, 229)
SEED_GROUPS = tuple(SEEDS[offset:offset + 5] for offset in range(0, 25, 5))


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
        b = baseline[indices][:, groups] if baseline.ndim == 2 else baseline[indices]
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
            endpoint_root = LOW2D_ENDPOINT if pair_spec["hard_low_2d"] else NONLOW2D_ENDPOINT
            path = endpoint_root / pair / f"group_{group}" / "result.json"
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
        lsalign_groups = []
        lsalign_pb = []
        lsalign_joint = []
        for group_index, seeds in enumerate(SEED_GROUPS, start=1):
            candidates = []
            for seed in seeds:
                meta_path = LSALIGN_GENERATION / pair / f"seed_{seed}" / "metadata.json"
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if (meta.get("status") != "success" or meta.get("pair") != pair or
                        meta.get("seed") != seed or meta.get("query_native_ligand_read") is not False):
                    raise RuntimeError("LS-align generation identity mismatch: " + str(meta_path))
                metric = meta["selected_metrics"]
                key_select = (-min(metric["pc_score_query"], metric["pc_score_template"]),
                              -metric["jaccard"], metric["ls_rmsd_A"], metric["mode"], seed)
                candidates.append((key_select, seed))
                input_hashes.append(sha256(meta_path))
            selected_seed = min(candidates)[1]
            endpoint_path = LSALIGN_ENDPOINT / pair / f"seed_{selected_seed}" / "result.json"
            endpoint = json.loads(endpoint_path.read_text(encoding="utf-8"))
            if (endpoint.get("status") != "success" or endpoint.get("pair") != pair or
                    endpoint.get("seed") != selected_seed or endpoint.get("generator_reads_native") is not False):
                raise RuntimeError("LS-align endpoint identity mismatch: " + str(endpoint_path))
            lsalign_groups.append(float(endpoint["rmsd_A"]))
            pb_path = LSALIGN_PB / pair / f"group_{group_index}" / "result.json"
            pb_result = json.loads(pb_path.read_text(encoding="utf-8"))
            if (pb_result.get("status") != "success" or pb_result.get("pair") != pair or
                    pb_result.get("group_index") != group_index or
                    pb_result.get("selected_seed") != selected_seed or
                    pb_result.get("selector_reads_native") is not False):
                raise RuntimeError("LS-align PB identity mismatch: " + str(pb_path))
            lsalign_pb.append(bool(pb_result["PB20_pass"]))
            lsalign_joint.append(bool(pb_result["joint_success"]))
            input_hashes.extend((sha256(endpoint_path), sha256(pb_path)))
        roshambo_groups = []
        roshambo_pb = []
        roshambo_joint = []
        for seed in (101, 127, 151, 179, 199):
            endpoint_path = ROSHAMBO_ENDPOINT / pair / f"seed_{seed}" / "result.json"
            endpoint = json.loads(endpoint_path.read_text(encoding="utf-8"))
            if (endpoint.get("status") != "success" or endpoint.get("pair") != pair or
                    endpoint.get("seed") != seed or endpoint.get("generator_reads_native") is not False):
                raise RuntimeError("ROSHAMBO endpoint identity mismatch: " + str(endpoint_path))
            pb_path = ROSHAMBO_PB / pair / f"seed_{seed}" / "result.json"
            pb_result = json.loads(pb_path.read_text(encoding="utf-8"))
            if (pb_result.get("status") != "success" or pb_result.get("pair") != pair or
                    pb_result.get("seed") != seed or
                    pb_result.get("generator_reads_native") is not False or
                    pb_result.get("target_native_ligand_read_for_physical_checks") is not False):
                raise RuntimeError("ROSHAMBO PB identity mismatch: " + str(pb_path))
            roshambo_groups.append(float(endpoint["rmsd_A"]))
            roshambo_pb.append(bool(pb_result["PB20_pass"]))
            roshambo_joint.append(bool(pb_result["joint_success"]))
            input_hashes.extend((sha256(endpoint_path), sha256(pb_path)))
        records.append({**row, **pair_spec, "groups": values, "pb": pb, "joint": joint,
                        "baselines": {"Iktos-3DAlign": float(align[key]), "Iktos-3DDock": float(dock[key]),
                                      "TEMPL-author": float(templ["rmsd_A"]),
                                      "LS-align-official-rerun": lsalign_groups,
                                      "ROSHAMBO-JCIM2024-rerun": roshambo_groups},
                        "lsalign_pb": lsalign_pb, "lsalign_joint": lsalign_joint,
                        "roshambo_pb": roshambo_pb, "roshambo_joint": roshambo_joint})
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
    result = {"analysis_version": "hybrid-conditional-v9", "input_digest": digest,
              "protocol": {"method_replicates": 5, "replicate_unit": "disjoint five-seed groups",
                           "success_threshold_A": 2.0, "difficulty_strata_overlap": True,
                           "prediction_selector_reads_native": False,
                           "hybrid_route": "hard_low_2d uses locked trust125 medoid; other rows use PB-valid receptor-score/consensus v2",
                           "v2_development_scope": "consensus weight selected using seed group 1; groups 2-5 are confirmation",
                           "baseline_runs_per_pair": {"Iktos-3DAlign": 1, "Iktos-3DDock": 1,
                                                      "TEMPL-author": 1, "LS-align-official-rerun": 5,
                                                      "ROSHAMBO-JCIM2024-rerun": 5},
                           "timeout_fallback": "2ica__3m6f used a preregistered deterministic non-random embedding fallback after the full schedule exceeded 900 seconds",
                           "baseline_sources": {"Iktos-3DAlign": "official same-row public output, commit 747d842",
                                                "Iktos-3DDock": "official same-row public output, commit 747d842",
                                                "TEMPL-author": "self-rerun author code, commit 727fd30, same rows",
                                                "LS-align-official-rerun": "self-rerun official executable, five disjoint five-seed groups, same rows",
                                                "ROSHAMBO-JCIM2024-rerun": "self-rerun JCIM 2024 author code, five seeds x ten conformers, same rows"},
                           "venue_restricted_exact_baselines": ["TEMPL-author", "LS-align-official-rerun", "ROSHAMBO-JCIM2024-rerun"],
                           "auxiliary_nonvenue_exact_baselines": ["Iktos-3DAlign", "Iktos-3DDock"]},
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
        eligible_dominance = []
        for method in ("Iktos-3DAlign", "Iktos-3DDock", "TEMPL-author", "LS-align-official-rerun", "ROSHAMBO-JCIM2024-rerun"):
            baseline = np.asarray([row["baselines"][method] for row in rows], dtype=float)
            baseline_pair_mean = baseline.mean(axis=1) if baseline.ndim == 2 else baseline
            baseline_group_rmsd = baseline.mean(axis=0) if baseline.ndim == 2 else np.repeat(baseline.mean(), 5)
            baseline_group_success = ((baseline <= 2.0).mean(axis=0) if baseline.ndim == 2
                                      else np.repeat((baseline <= 2.0).mean(), 5))
            pair_delta = proposed.mean(axis=1) - baseline_pair_mean
            signflip = exact_family_signflip(rows, pair_delta)
            raw_p.append((name, method, signflip["one_sided_p"]))
            all_groups = [bool(group_success[i] >= baseline_group_success[i] and
                               group_rmsd[i] < baseline_group_rmsd[i])
                          for i in range(5)]
            baseline_pb = (np.asarray([row["lsalign_pb"] for row in rows], dtype=bool)
                           if method == "LS-align-official-rerun" else
                           np.asarray([row["roshambo_pb"] for row in rows], dtype=bool)
                           if method == "ROSHAMBO-JCIM2024-rerun" else None)
            baseline_joint = (np.asarray([row["lsalign_joint"] for row in rows], dtype=bool)
                              if method == "LS-align-official-rerun" else
                              np.asarray([row["roshambo_joint"] for row in rows], dtype=bool)
                              if method == "ROSHAMBO-JCIM2024-rerun" else None)
            entry["baselines"][method] = {
                "success": float((baseline <= 2.0).mean()), "success_std": sample_std(baseline_group_success) if baseline.ndim == 2 else None,
                "PB20": float(baseline_pb.mean()) if baseline_pb is not None else None,
                "PB20_std": sample_std(baseline_pb.mean(axis=0)) if baseline_pb is not None else None,
                "joint_success": float(baseline_joint.mean()) if baseline_joint is not None else None,
                "joint_success_std": sample_std(baseline_joint.mean(axis=0)) if baseline_joint is not None else None,
                "mean_rmsd_A": float(baseline.mean()), "rmsd_A_std": sample_std(baseline_group_rmsd) if baseline.ndim == 2 else None,
                "per_group_mean_rmsd_A": baseline_group_rmsd.tolist(),
                "delta_proposed_mean_rmsd_A": float(proposed.mean() - baseline.mean()),
                "pair_win_rate": float((proposed.mean(axis=1) < baseline_pair_mean).mean()),
                "all_five_groups_dominate": all(all_groups), "per_group_dominance": all_groups,
                "all_25_cross_replicates_dominate": bool(
                    group_rmsd.max() < baseline_group_rmsd.min()
                    and group_success.min() >= baseline_group_success.max()),
                "hierarchical_bootstrap": bootstrap(rows, proposed, baseline),
                "family_signflip": signflip,
            }
            dominance.append(all(all_groups))
            if method in ("TEMPL-author", "LS-align-official-rerun", "ROSHAMBO-JCIM2024-rerun"):
                eligible_dominance.append(all(all_groups))
        entry["same_row_all_five_group_dominance"] = all(dominance)
        entry["venue_restricted_all_five_group_dominance"] = all(eligible_dominance)
        entry["claim_scope"] = ("repeated venue-restricted same-row leader"
                                if all(eligible_dominance) else "no venue-restricted dominance claim")
        result["strata"][name] = entry
    adjusted = holm_adjust(raw_p)
    for name, method, _ in raw_p:
        result["strata"][name]["baselines"][method]["family_signflip"]["holm_adjusted_p_across_35"] = adjusted[(name, method)]
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
        with (stage / "main_table.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["stratum", "n_pairs", "n_families", "method", "mean_rmsd_A", "rmsd_std_across_5_groups", "success", "PB20", "joint_success", "source"])
            for name, entry in result["strata"].items():
                proposed = entry["proposed"]
                writer.writerow([name, entry["n_pairs"], entry["n_families"], "RouteAlign", proposed["rmsd_A_mean"], proposed["rmsd_A_std"], proposed["success_mean"], proposed["PB20_mean"], proposed["joint_mean"], "self-rerun"])
                for method, baseline in entry["baselines"].items():
                    if method.startswith("Iktos"):
                        source = "official public technical-report output; auxiliary"
                    elif method.startswith("LS-align"):
                        source = "self-rerun official executable; five seed groups"
                    elif method.startswith("ROSHAMBO"):
                        source = "self-rerun JCIM 2024 author code; five seeds x ten conformers"
                    else:
                        source = "self-rerun author code"
                    writer.writerow([name, entry["n_pairs"], entry["n_families"], method, baseline["mean_rmsd_A"], baseline.get("rmsd_A_std") or "", baseline["success"], baseline.get("PB20") if baseline.get("PB20") is not None else "", baseline.get("joint_success") if baseline.get("joint_success") is not None else "", source])
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
