# Current verification scope

The correction recalculates chemical graph distances for all 220 groups and
re-evaluates 660 selected endpoints from 1,100 unchanged candidate poses.
All 1,100 physical checks are rerun and agree with the archived pass flags.

Checks cover atom renumbering, aromatic and resonance symmetry, nonzero
receptor-frame translation, bounded timeout, failure cleanup and hash-sensitive
recovery. Observed candidate comparisons include up to 48 complete atom maps.

`python -m routealign.reproduce` verifies the corrected seven-stratum primary
analysis and three-stratum same-candidate selector ablation. Statistical
intervals use 20,000 hierarchical resamples; family tests and Holm adjustments
are recomputed. Independent external-model integrity review is unavailable in
this correction run; execution checks are not presented as such a review.

The additional input-stereochemistry audit remains an unresolved molecular
correctness finding, explicitly recorded with all affected pair and seed IDs.
