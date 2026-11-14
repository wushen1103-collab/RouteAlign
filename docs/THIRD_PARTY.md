# Data and software provenance

RouteAlign uses external scientific tools without bundling their executables,
training weights or full data archives. Downloaded materials retain their own
licences and notices.

| Component | Source / frozen version | Role |
|---|---|---|
| AlignDockBench | https://doi.org/10.5281/zenodo.15395813 | Original template/query structures; archive MD5 `ca8d9076efa4b4661c74c47decee2ff1`; CC BY-NC-ND 4.0 |
| RDKit | https://www.rdkit.org/ ; 2026.09.1 | Main conformer generation, graph operations and endpoint chemistry |
| TEMPL-generation RDKit | 2024.9.6 with NumPy 1.26.4 | Separate environment matching the original TEMPL rerun |
| PoseBusters | https://github.com/maabuu/posebusters ; 0.5.1 | Physical checks and symmetry-aware RMSD |
| Smina | https://sourceforge.net/projects/smina/ | Score-only evaluation and validity-triggered minimization |
| TEMPL | https://github.com/fulopjoz/templ-pipeline ; `727fd306e4066aa39b05a6e51e7833c86d90df98` | Rerun author-code baseline; MIT source notice downloaded with the code |
| LS-align | https://zhanggroup.org/LS-align/ | Rerun author executable, rigid and flexible alignment modes |
| Open Babel | https://openbabel.org/ ; 3.1.1 | LS-align SDF/MOL2 conversion |
| ROSHAMBO | https://github.com/molecularinformatics/roshambo ; `916fe97f81a3fdb48e7bb85384eb597836d2b01e` | Rerun author code in a separate RDKit 2023.09.6 environment; upstream GPL-3.0 |
| Iktos engine | https://github.com/iktos/engine ; `747d842` | Two auxiliary published-output CSVs; MIT |

The original Smina executable SHA-256 is
`ffe5e1e78c947f76b0df8805e2c54383d0bbaf2e827a633b643a708cf682a958`.
The LS-align executable SHA-256 is
`598afc0cd5f8994d9aba159b45315de12fa5ecdff4a226de07b281959c9bf1bc`.
These hashes identify software inputs, not user machines.

The TEMPL wrapper adapts its two coordinate-map assignments to RDKit's
`SetCoordMap` API. The ROSHAMBO wrapper caps upstream multiprocessing to the
explicit CPU allocation and restores poses to the template coordinate frame.
Those compatibility changes do not introduce native query coordinates into
generation or selection.

Scientific method references are available through the linked source projects
and their publications. No third-party authorship is represented as RouteAlign
authorship, and the RouteAlign licence does not relicense third-party materials.
