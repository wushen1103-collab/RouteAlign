# RouteAlign

Input-conditioned consensus for template-based molecular pose prediction.

RouteAlign starts from a query molecular graph, a template ligand pose and its
receptor. Common-core grafting generates seed-level candidates. Physical-validity
checks trigger a bounded repair step. A fixed Morgan similarity threshold routes
low-similarity inputs to candidate-medoid selection and other inputs to
receptor-score/consensus selection. The selected coordinates are saved before
native-pose RMSD is evaluated.

## Reproduce the reported numbers

The small frozen observation archive includes all 44 pairs, 25 template families,
five disjoint groups of five seeds, three rerun baseline families and two
separately labelled published-output baselines. It includes a three-rule, identical-candidate selector ablation. No GPU or dataset download is needed to recompute the tables and tests.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e . -i https://pypi.tuna.tsinghua.edu.cn/simple
python -m routealign.reproduce
```

The command recomputes the corrected primary analysis, all seven overlapping strata, and
the same-candidate selector ablation, with hierarchical
bootstrap intervals, family-level sign-flip tests and Holm corrections. It checks
the regenerated numerical results against the frozen reference summaries and
writes machine-readable JSON and CSV tables to `outputs/statistics/`. This is a
statistics-level reproduction; the molecular pipeline below regenerates poses.

For the 44-pair primary cohort, the frozen RouteAlign result is mean RMSD
0.6877 ± 0.0176 Å, RMSD ≤ 2 Å success 99.55 ± 1.02%, physical-validity pass
97.73 ± 0.00%, and joint success 97.27 ± 1.02%. The standard deviations are across
five independent seed groups, not individual targets. TEMPL has one author-code
run per pair; its uncertainty is not represented as five independent runs.

## Correspondence correction and selector evidence

Candidate distances now minimize over full chemical graph mappings, including
graph automorphisms and resonance-equivalent atoms, without fitting coordinates.
On identical candidates, mean RMSD is 0.693556 A for uniform medoid, 0.687802 A
for uniform receptor-consensus, and 0.687659 A for routing. The routed-minus-
receptor difference has a 95% hierarchical-bootstrap interval of
[-0.000621, 0.000199] A and Holm-adjusted p = 1.000 across six ablation tests.
These controls do not establish an accuracy advantage from routing.

An input-graph audit additionally records 40 candidate and 10 selected-output
stereochemical-annotation discrepancies across three pairs. They are retained
in `reference/stereochemistry_audit.json`; the 19-check physical endpoint does
not replace this audit. Molecular correctness requires this separate review.

Run `python -m routealign.ablation` for the selector-only statistical replay.
The older `uniform_summary.json` describes a historical whole-system comparison,
which is no longer used as evidence isolating final routing.

## Regenerate molecular predictions

The molecular workers run on Linux. Install the chemistry extras in the active
Python 3.10 environment, then prepare a workspace using the original dataset:

```bash
python -m pip install -e '.[chemistry]' -i https://pypi.tuna.tsinghua.edu.cn/simple
python -m routealign.prepare --download
python -m routealign.tools --smina-download
```

An existing extracted dataset can be reused with
`python -m routealign.prepare --dataset /path/to/AlignDockBench`.
An existing matching Smina binary can be reused with
`python -m routealign.tools --smina /path/to/smina`.
Downloads and binaries are excluded from Git.

Start with one task and inspect its outputs:

```bash
python -m routealign.run --stage graft --index 300
python -m routealign.run --stage repair --index 300
python -m routealign.run --stage bounded_repair --index 300
```

Continue with the full stage sequence in [docs/REPRODUCING.md](docs/REPRODUCING.md).
Every task has a finite timeout, an isolated process group, atomic final outputs
and input-sensitive resume checks. Defaults are one worker, one CPU per task and
hidden GPUs. On machines with more than 30 allowed CPUs, the runner retains at
least 30 CPUs when accepting a larger `--workers` value. ROSHAMBO generation
requires a separate environment, an explicit GPU and one worker.

## Contents

| Directory | Contents |
|---|---|
| `routealign/` | Portable scientific workers, resource-controlled entry points and statistics |
| `config/` | Exact cohort identities, input descriptors, seed groups and frozen protocol |
| `reference/` | Compressed observations, numerical summaries and checksums |
| `docs/` | Full reproduction sequence, endpoint definitions and third-party provenance |
| `tests/` | CPU/device isolation and process-group timeout tests |

The supporting CSVs include the endpoint-level observations, physical-check
breakdown, continuous geometry, repair changes, atom-order stress checks and
structural-case values. Their fields and aggregation units are described in
[docs/SUPPORTING_DATA.md](docs/SUPPORTING_DATA.md).

## Data and third-party methods

The original AlignDockBench archive is available at
[Zenodo, DOI 10.5281/zenodo.15395813](https://doi.org/10.5281/zenodo.15395813).
It is downloaded separately and retains its original CC BY-NC-ND 4.0 licence.
See [docs/THIRD_PARTY.md](docs/THIRD_PARTY.md) for software sources and the
distinction between rerun and published-output comparisons.

The MIT licence covers the RouteAlign code. It does not replace the licences of
the original dataset, external executables or downloaded baseline source code.
