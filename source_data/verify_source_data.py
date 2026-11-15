"""Verify the versioned Source Data with one CPU and Python's standard library."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys


def boolean(value):
    if value in [True, "True", "true", "1"]:
        return True
    if value in [False, "False", "false", "0"]:
        return False
    raise ValueError("Invalid Boolean value")


def verify(directory):
    directory = directory.resolve()
    manifest = json.loads((directory / "checksums.json").read_text(encoding="utf-8"))
    for name, info in manifest["files"].items():
        path = (directory / name).resolve()
        if not path.is_relative_to(directory) or not path.is_file():
            raise ValueError("Missing or unsafe payload entry: " + name)
        data = path.read_bytes()
        if len(data) != info["bytes"] or hashlib.sha256(data).hexdigest() != info["sha256"]:
            raise ValueError("Checksum mismatch: " + name)

    def rows(name):
        with (directory / name).open(encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))

    checks = 0

    def check(condition, message):
        nonlocal checks
        checks += 1
        if not condition:
            raise ValueError(message)

    def close(a, b, message):
        check(math.isclose(float(a), float(b), rel_tol=1e-11, abs_tol=1e-12), message)

    cohort = rows("cohort.csv")
    check(len(cohort) == 44 and len({r["family"] for r in cohort}) == 25, "Cohort size")
    check(len(rows("eligibility.csv")) == 369, "Eligibility population")
    endpoints = rows("endpoints.csv")
    check(len(endpoints) == 704, "Endpoint count")
    route = {(r["pair"], int(r["repeat_group"])): r for r in endpoints if r["method"] == "RouteAlign"}
    check(len(route) == 220, "RouteAlign endpoint count")
    selected = {(r["pair"], int(r["repeat_group"]), r["rule"]): r for r in rows("selector_endpoints.csv")}
    check(len(selected) == 660, "Selector endpoint count")
    groups = json.loads(gzip.decompress((directory / "candidate_groups.json.gz").read_bytes()))
    check(len(groups) == 220, "Five-candidate group count")

    def zscore(values):
        sd = statistics.pstdev(values)
        return [0.0] * len(values) if sd <= 1e-12 else [(v - statistics.fmean(values)) / sd for v in values]

    units = []
    for g in groups:
        f = g["features"]
        d = f["distance_matrix"]
        check(len(d) == 5 and all(len(row) == 5 for row in d), "Distance-matrix shape")
        for i in range(5):
            close(d[i][i], 0, "Distance-matrix diagonal")
            for j in range(i):
                close(d[i][j], d[j][i], "Distance-matrix symmetry")
        med = [statistics.median([d[i][j] for j in range(5) if j != i]) for i in range(5)]
        avg = [statistics.fmean([d[i][j] for j in range(5) if j != i]) for i in range(5)]
        seeds = f["seeds"]
        mi = min(range(5), key=lambda i: (med[i], avg[i], seeds[i]))
        zs, zm = zscore(f["receptor_scores"]), zscore(med)
        eligible = [i for i in range(5) if f["physical_pass"][i]] or list(range(5))
        ri = min(eligible, key=lambda i: (zs[i] + 5 * zm[i], seeds[i]))
        choices = {"medoid": mi, "receptor_consensus": ri, "routed": mi if f["low_similarity"] else ri}
        for rule, i in choices.items():
            row = selected[(g["pair"], g["group"], rule)]
            check(i == g["choices"][rule]["index"] == int(row["candidate_index"]), "Selector choice")
            check(int(row["seed"]) == seeds[i], "Selected seed")
            check(row["pose_sha256"] == g["candidate_hashes"][i], "Selected structure identity")
            rmsd = g["endpoints"]["rmsd_A"][i]
            pb = f["physical_pass"][i]
            close(row["rmsd_A"], rmsd, "Selected endpoint RMSD")
            check(boolean(row["success_2A"]) == (rmsd <= 2), "Accuracy threshold")
            check(boolean(row["physical_pass"]) == pb, "Physical endpoint")
            check(boolean(row["joint_success"]) == (rmsd <= 2 and pb), "Joint endpoint")
            units.append({"group": g["group"], "low": f["low_similarity"], "rule": rule, "rmsd_A": rmsd, "success_pct": 100 * (rmsd <= 2), "physical_pass_pct": 100 * pb, "joint_success_pct": 100 * (rmsd <= 2 and pb)})
            if rule == "routed":
                row = route[(g["pair"], g["group"])]
                close(row["rmsd_A"], rmsd, "RouteAlign export agrees")
                check(boolean(row["physical_pass"]) == pb, "RouteAlign physical export agrees")

    for row in rows("table_s8_ablation.csv"):
        subset = [r for r in units if r["rule"] == row["rule"] and (row["stratum"] == "all" or r["low"] == (row["stratum"] == "low"))]
        for metric in ["rmsd_A", "success_pct", "physical_pass_pct", "joint_success_pct"]:
            values = [statistics.fmean(r[metric] for r in subset if r["group"] == g) for g in range(1, 6)]
            close(row[metric + "_mean"], statistics.fmean(values), "Ablation mean")
            close(row[metric + "_sd"], statistics.stdev(values), "Ablation sample SD")

    for row in rows("table_1.csv"):
        subset = [r for r in endpoints if r["method"] == row["method"]]
        for metric, source in [("rmsd_A", "rmsd_A"), ("success_pct", "success_2A"), ("physical_pass_pct", "physical_pass"), ("joint_success_pct", "joint_success")]:
            if row[metric + "_mean"] == "":
                continue
            count = int(row["replicate_groups"])
            values = [statistics.fmean(float(r[source]) if metric == "rmsd_A" else 100 * boolean(r[source]) for r in subset if int(r["repeat_group"]) == g) for g in range(1, count + 1)]
            close(row[metric + "_mean"], statistics.fmean(values), "Primary table mean")
            if count == 5:
                close(row[metric + "_sd"], statistics.stdev(values), "Primary table sample SD")
            else:
                check(row[metric + "_sd"] == "", "Single-output SD remains empty")

    means = {(pair, method): statistics.fmean(float(r["rmsd_A"]) for r in endpoints if r["pair"] == pair and r["method"] == method) for pair in {r["pair"] for r in endpoints} for method in ["RouteAlign", "TEMPL", "LS-align", "ROSHAMBO"]}
    for row in rows("figure_2.csv"):
        close(row["RouteAlign_pair_mean_rmsd_A"], means[(row["pair"], "RouteAlign")], "Fig. 2 proposed value")
        close(row["baseline_pair_mean_rmsd_A"], means[(row["pair"], row["baseline"])], "Fig. 2 comparator value")
    for row in rows("figure_s2.csv"):
        close(row["pair_mean_rmsd_A"], means[(row["pair"], row["method"])], "Fig. S2 value")
    for row in rows("figure_s5.csv"):
        close(row["RouteAlign_rmsd_A"], route[(row["pair"], int(row["group"]))]["rmsd_A"], "Fig. S5 value")
    for row in rows("physical_checks.csv"):
        names = [k for k in row if k not in ["pair", "family", "group", "method", "physical_pass"]]
        check(len(names) == 19, "Physical-check count")
        check(all(boolean(row[k]) for k in names) == boolean(row["physical_pass"]), "Physical-check conjunction")

    current = list(route.values())
    close(statistics.fmean(float(r["rmsd_A"]) for r in current), 0.6729458703064845, "Corrected headline RMSD")
    check(sum(boolean(r["success_2A"]) for r in current) == 219, "219 accuracy successes")
    check(sum(boolean(r["physical_pass"]) for r in current) == 215, "215 physical successes")
    check(sum(boolean(r["joint_success"]) for r in current) == 214, "214 joint successes")
    return {"status": "PASS", "version": manifest["version"], "payload_files": len(manifest["files"]), "scientific_checks": checks, "pairs": 44, "replicate_groups": 5, "candidate_choices": 660, "mean_rmsd_A": statistics.fmean(float(r["rmsd_A"]) for r in current)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    here = Path(__file__).resolve()
    default = here.parent if here.name == "verify_source_data.py" else here.parents[1] / "source_data"
    parser.add_argument("--directory", type=Path, default=default)
    args = parser.parse_args()
    try:
        result = verify(args.directory)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(json.dumps({"status": "FAIL", "reason": str(exc)}), file=sys.stderr)
        raise SystemExit(1)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
