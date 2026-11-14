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
separately labelled published-output baselines. It includes the uniform-medoid
ablation. No GPU or dataset download is needed to recompute the tables and tests.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e . -i https://pypi.tuna.tsinghua.edu.cn/simple
python -m routealign.reproduce
```

The command recomputes both analyses, all seven overlapping strata, hierarchical
bootstrap intervals, family-level sign-flip tests and Holm corrections. It checks
the regenerated numerical results against the frozen reference summaries and
writes machine-readable JSON and CSV tables to `outputs/statistics/`. This is a
statistics-level reproduction; the molecular pipeline below regenerates poses.

For the 44-pair primary cohort, the frozen RouteAlign result is mean RMSD
0.6869 ± 0.0185 Å, RMSD ≤ 2 Å success 99.55 ± 1.02%, physical-validity pass
97.73 ± 0.00%, and joint success 97.27 ± 1.02%. The standard deviations are across
five independent seed groups, not individual targets. TEMPL has one author-code
run per pair; its uncertainty is not represented as five independent runs.

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
