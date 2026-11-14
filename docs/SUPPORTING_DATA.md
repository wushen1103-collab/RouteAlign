# Supporting numerical data

All files below are in `reference/supporting/`. They contain scientific
identifiers and measured values, without execution paths or machine metadata.

| File | Rows | Meaning |
|---|---:|---|
| `endpoint_observations.csv` | 704 | Per-pair/group RouteAlign, LS-align and ROSHAMBO endpoints, plus one TEMPL endpoint per pair |
| `family_effects.csv` | 21 | Three rerun comparators × seven overlapping strata; RMSD difference, interval, wins and corrected test |
| `physical_checks.csv` | 57 | Nineteen PoseBusters checks × three five-group methods; failures and denominators |
| `continuous_geometry.csv` | 2,640 | Selected-pose geometry diagnostics, retaining feature names, display scale and missing-value flags |
| `repair_observations.csv` | 110 | Selected observations for which repair was triggered, with before/after diagnostic values |
| `atom_order_stress.csv` | 9 | Distinct atom-element sequences and mixed five-candidate groups in the ordering stress check |
| `structural_cases.csv` | 4 | Case identity, selected seed/group, RMSD, physical pass, repair movement and coordinate checksums |

`physical_checks.csv` gives counts before percentage formatting. Each method has
220 selected outputs (44 pairs × five groups); the 19 checks are not independent
experimental replicates. The named features in `continuous_geometry.csv` are
coordinate diagnostics; they do not replace the full physical-validity endpoint.

Repair rows describe the subset that triggered a repair attempt. An accepted
repair must be interpreted together with its acceptance flag and movement; these
rows are not a separately sampled benchmark. `structural_cases.csv` identifies
the native and predicted coordinate files by content checksum without bundling
the original dataset structures.

`routealign.reproduce` regenerates the primary and uniform-ablation summaries,
including intervals and tests, directly from the frozen per-pair observations.
The supporting geometry CSVs are the archived deterministic diagnostics from
those saved coordinates. Dataset downloads, fresh pose generation, and frozen
observation aggregation are distinct reproduction operations.
