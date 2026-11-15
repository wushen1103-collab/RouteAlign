# Supporting numerical data

The canonical corrected figure/table inputs are now in
[`source_data/`](../source_data/README.md), version `source-data-v1.0.1`.
That snapshot supplies the explicit display-number mapping, units, regenerated
current endpoint tables and an independently downloadable archive. Verify it with
`python -m routealign.source_data` or the verifier bundled in the ZIP.
The `reference/supporting/` files below retain the statistical/diagnostic export
history used by existing reproduction commands; use the versioned Source Data
directory when locating values for the final figures and tables.

All files below are in `reference/supporting/`. They contain scientific
identifiers and measured values, without execution paths or machine metadata.

| File | Rows | Meaning |
|---|---:|---|
| `endpoint_observations.csv` | 1,409 | Per-pair/group selector, RouteAlign, LS-align and ROSHAMBO endpoints, plus one TEMPL endpoint per pair |
| `family_effects.csv` | 43 | Comparator and selector effects over the prespecified strata; RMSD difference, interval, wins and corrected test |
| `physical_checks.csv` | 57 | Nineteen PoseBusters checks × three five-group methods; failures and denominators |
| `continuous_geometry.csv` | 2,640 | Selected-pose geometry diagnostics, retaining feature names, display scale and missing-value flags |
| `repair_observations.csv` | 106 | Selected observations for which repair was triggered, with before/after diagnostic values |
| `atom_order_stress.csv` | 9 | Distinct atom-element sequences and mixed five-candidate groups in the ordering stress check |
| `structural_cases.csv` | 4 | Case identity, selected seed/group, RMSD, physical pass, repair movement and coordinate checksums |
| `same_candidate_ablation.csv` | 90 | Three selector rules on the same corrected candidate pool, split by similarity regime |

`physical_checks.csv` gives counts before percentage formatting. Each method has
220 selected outputs (44 pairs × five groups); the 19 checks are not independent
experimental replicates. The named features in `continuous_geometry.csv` are
coordinate diagnostics; they do not replace the full physical-validity endpoint.

Repair rows describe the subset that triggered a repair attempt. An accepted
repair must be interpreted together with its acceptance flag and movement; these
rows are not a separately sampled benchmark. `structural_cases.csv` identifies
the native and predicted coordinate files by content checksum without bundling
the original dataset structures.

`routealign.reproduce` regenerates the corrected primary and identical-candidate ablation summaries,
including intervals and tests, directly from the frozen per-pair observations.
The supporting geometry CSVs are the archived deterministic diagnostics from
those saved coordinates. Dataset downloads, fresh pose generation, and frozen
observation aggregation are distinct reproduction operations.

The active selector comparison uses `same_candidate_groups.json.gz`, with
`ablation_summary.json` and `supporting/same_candidate_ablation.csv`. Endpoint
observations contain all three rules (660 selected endpoints). Historical
uniform-system observations are retained and are not final-routing controls.
