# Reproduction guide

## 1. Environment and inputs

Use Python 3.10 with NumPy 1.26.4, RDKit 2026.09.1 and PoseBusters 0.5.1 for the
molecular pipeline. Install `.[chemistry]` from the repository root. Frozen
statistics require only NumPy and also run on newer Python versions.

`routealign.prepare` copies the required template ligands/receptors and native
query ligands from the original AlignDockBench archive into the workspace. It
verifies every ligand against the frozen cohort's SHA-256. Generation and
selection read query graphs from `config/cohort.csv`; native query coordinates
are read by endpoint evaluation after selected coordinates have been saved.

Choose a separate working directory by passing `--workspace /path/to/workspace`
to preparation, tools and run commands. Set `ROUTEALIGN_WORKSPACE` to that same
directory when using `routealign.reproduce --live`.

## 2. Verify control before scaling up

```bash
python -m unittest discover -s tests -v
python -m routealign.run --stage graft --index 300 --fault-after-candidates
```

The injected failure returns 97 without committing a task output. A clean run
computes the task; a repeated identical run resumes without replacing it:

```bash
python -m routealign.run --stage graft --index 300
python -m routealign.run --stage graft --index 300
```

Use `--timeout SECONDS` to impose a deadline. The runner terminates the process
group, including child tools. Fault and hold switches are available on workers
at their commit boundaries. A changed input is rejected instead of silently
reusing an incompatible output.

Default concurrency is one task on one CPU. `--workers N` is accepted only within
the allowed affinity mask and CPU reserve. GPU generation is restricted to a
single explicitly requested ROSHAMBO task at a time; its worker also checks
device identity/load and takes an exclusive per-device lock. An idle-device
reading does not start additional jobs.

## 3. Regenerate RouteAlign and the historical uniform system

There are 44 pairs × 25 seeds = 1,100 seed-level tasks. Seeds form five disjoint
groups of five. The full protocol and pair order are in `config/protocol.json`.
The deterministic fallback for `2ica_2IC_A_1__3m6f` is selected automatically by
the runner; it retains the original non-random embedding attempt schedule.

```bash
python -m routealign.run --stage graft --all
python -m routealign.run --stage repair --all
python -m routealign.run --stage bounded_repair --all
python -m routealign.run --stage medoid_select --all
python -m routealign.run --stage medoid_evaluate --all
```

These stages provide the low-similarity route and the uniform-medoid ablation.
The receptor-conditioned route covers the other 18 pairs, giving 450 seed tasks
and 90 group tasks:

```bash
python -m routealign.run --stage receptor_score --all
python -m routealign.run --stage receptor_repair --all
python -m routealign.run --stage receptor_select --all
python -m routealign.run --stage receptor_evaluate --all
```

Tasks within a stage can be restarted independently. Complete each stage before
starting a dependent stage. For a subset, use `--index I` or `--indices I,J,K`.
Graft/repair indices are `pair_index * 25 + seed_position`. Medoid group indices
are `pair_index * 5 + group_position`, with zero-based positions. Receptor-route
indices use the protocol order filtered to the 18 non-low-similarity pairs.

## 4. Rerun the matched baselines

### TEMPL

Use a separate Python 3.10 environment with RDKit 2024.9.6 and NumPy 1.26.4
for TEMPL generation (`pip install -e '.[templ]'`). The reported TEMPL generator
and common endpoint evaluator deliberately use different RDKit versions.

```bash
# In the TEMPL environment:
python -m routealign.tools --templ
python -m routealign.run --stage templ_generate --all
# Switch back to the common endpoint environment:
python -m routealign.run --stage templ_evaluate --all
```

The wrapper downloads the frozen author source snapshot, checks its checksums
and uses the same coordinate-map compatibility change as the reported run.
It produces one output per pair. It does not manufacture five independent runs.

### LS-align

Obtain the author executable from https://zhanggroup.org/LS-align/ and install
Open Babel 3.1.1 in the active environment. The wrapper records the actual Open
Babel binary checksum; its installation path is configurable through
`ROUTEALIGN_OBABEL`. The LS-align executable is checked against the reported
SHA-256 when copied into the workspace.

```bash
python -m routealign.tools --lsalign /path/to/LSalign
python -m routealign.run --stage lsalign_generate --all
python -m routealign.run --stage lsalign_evaluate --all
python -m routealign.run --stage lsalign_validate --all
```

This evaluates the same 25 seeds per pair. The selected member of each five-seed
group is chosen from LS-align scores before native-pose evaluation.

### ROSHAMBO

ROSHAMBO uses a separate environment with RDKit 2023.09.6 and the author source
snapshot recorded in [THIRD_PARTY.md](THIRD_PARTY.md). Build/install the author
CUDA extension in that environment. Put the checkout at
`WORKSPACE/third_party/roshambo`. Its source commit and setup-file checksum are
checked; the locally built extension checksum is recorded without requiring a
binary compiled on another machine.

```bash
# Run generation in the ROSHAMBO environment. Replace GPU_ID explicitly.
python -m routealign.run --stage roshambo_generate --all --gpu GPU_ID
# Switch back to the common endpoint environment for these stages.
python -m routealign.run --stage roshambo_evaluate --all
python -m routealign.run --stage roshambo_validate --all
```

There are 220 tasks: five seeds per pair, ten conformers requested per seed.
The reference-frame restoration and upstream worker-pool cap are included in
the wrapper. See `config/roshambo_protocol.md` for the exact settings.

### Published-output comparators

Iktos-3DAlign and Iktos-3DDock are auxiliary same-row published-output references.
They are labelled separately from self-rerun methods in every exported table.
Fetch the frozen public CSVs for live aggregation:

```bash
python -m routealign.tools --published-baselines
```

Their original MIT licence is downloaded alongside the CSVs. Full-benchmark
literature results from other protocols are not substituted for matched-cohort
reruns in the frozen observation archive.

## 5. Aggregate and compare

```bash
# Recompute the frozen experimental record, without generating poses.
python -m routealign.reproduce
# After rerunning all method and baseline stages, aggregate the new outputs.
python -m routealign.reproduce --live --output outputs/live_statistics
```

The live mode applies the same tests and checks the newly computed values
against the frozen reference. A mismatch raises an error with the affected
quantity instead of silently overwriting reference results. Different chemistry
or tool versions can alter generated conformers; use the recorded versions for
the reported experimental protocol.

## 6. Endpoints and statistics

RMSD is the symmetry-aware heavy-atom coordinate error relative to the native
query pose, in the supplied receptor frame. Success means RMSD ≤ 2 Å. Physical
pass means passing all checks returned by PoseBusters 0.5.1 with `config="dock"`;
joint success requires both physical pass and RMSD success. These definitions
are used consistently for the five-group summaries.

Mean ± sample SD is computed across five independent seed groups. Bootstrap
sampling resamples template families and seed groups, with 20,000 resamples.
Family-level sign flips use exact enumeration for at most 16 families and a
fixed-seed 200,000-draw Monte Carlo test otherwise. The 35 hybrid-analysis tests
use one Holm correction family. Difficulty strata overlap; their counts should
not be summed as independent datasets.

## Corrected selector-only control

Both final selectors use `routealign.atom_correspondence` for chemical graph
mapping, symmetry and receptor-frame distances. `python -m routealign.ablation`
replays all 220 groups and three rules on identical frozen candidates and
recomputes the paired statistics. This is the active final-selector ablation.
The historical uniform system also changed generation and cannot isolate routing.
SDF chirality is audited separately; coordinate matching preserves connectivity
and isotopes without changing or fitting coordinates. See the audit JSON.
