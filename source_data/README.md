# RouteAlign Source Data

Version **1.0** contains the corrected numerical evidence for 44 molecular pairs,
25 template families and five disjoint groups of five candidate seeds. It is
aligned with the graph- and stereochemistry-corrected workflow: mean RMSD
0.6729458703064845 A and joint success 97.27272727272729%.

## Download and verify

- Browse the versioned files: `https://github.com/wushen1103-collab/RouteAlign/tree/source-data-v1.0/source_data`
- Download the single archive: `https://raw.githubusercontent.com/wushen1103-collab/RouteAlign/source-data-v1.0/source_data/RouteAlign_Source_Data.zip`
- Verify after cloning: `python -m routealign.source_data`
- Verify an extracted archive: `python -m routealign.source_data --directory /path/to/extracted-data`
- Or run the verifier bundled in the extracted archive: `python verify_source_data.py`.
- Recompute full statistical tests from the repository reference observations:
  `python -m routealign.reproduce` and `python -m routealign.ablation`.

Verification of this source-data snapshot needs only the Python standard library,
one CPU and no GPU, dataset download, docking or geometry optimization.

## Files and display mapping

| Display or evidence | Data |
|---|---|
| Fig. 1 | Workflow schematic; molecular case identities in `figure_s6.csv` |
| Figs. 2-6 | `figure_2.csv` through `figure_6.csv` |
| Figs. S1-S7 | `figure_s1.csv` through `figure_s7.csv` |
| Table 1 | `table_1.csv` |
| Table S1 | `table_s1_methods.csv` |
| Table S2 | `cohort.csv` |
| Tables S3-S8 | `table_s3_strata.csv` through `table_s8_ablation.csv`, with descriptive suffixes |
| Algorithm S1 | `candidate_groups.json.gz`, `selector_endpoints.csv`, repository `config/protocol.json` |
| Full source-population eligibility | `eligibility.csv` |
| All matched-method observations | `endpoints.csv` |
| All 660 same-candidate selector outcomes | `selector_endpoints.csv` |
| All 19 physical checks per selected observation | `physical_checks.csv` |
| Correspondence and stereochemistry audits | `atom_order_audit.csv`, `graph_correspondence_changes.csv`, `stereochemistry_correction.json` |
| Units, columns and aggregation | `data_dictionary.json` |

CSV values retain numerical precision before manuscript rounding. Empty cells mean
unavailable or inapplicable. TEMPL has one output per pair, so its across-group SD
and unmeasured physical endpoint are empty. Percent columns are explicitly named.
The overlapping strata are not independent datasets. Geometry-feature observations
are not independent replicates of the 19-check endpoint.

## Coordinate-specific diagnostics

The continuous-geometry and triggered-repair files reproduce the frozen measurements
used in the displayed diagnostic panels. They retain their coordinate-specific
role. `diagnostic_pose_hashes.csv` identifies the checked coordinates and the final
selected candidate separately; a false byte-identity flag does not establish a
geometric difference. Use `endpoints.csv` and `selector_endpoints.csv` for the final
corrected accuracy/joint endpoints, and do not infer them from a geometry feature.
The 53 triggered physical-repair observations and the 40 later stereochemical
corrections are distinct operations and have separate records.

Structural-case measurements, common-core mappings and coordinate checksums are
included. The original input structures remain in AlignDockBench, Zenodo record
15395813 (DOI `10.5281/zenodo.15395813`), under its original CC BY-NC-ND 4.0 license.
They are not duplicated here. The repository molecular workflow regenerates poses.

## Provenance and licensing

`checksums.json` records every expanded payload file by SHA-256. The ZIP contains
that same payload and checksum file; the ZIP itself is not recursively hashed.
Snapshot manifests contain scientific identifiers only, with no execution host,
absolute paths, credentials or author/submission metadata.

This release's original numerical tables, documentation and RouteAlign software
follow the repository MIT License, copied to `LICENSE.txt`. Original benchmark
data and external tools retain their own licenses.
