# Current verification scope

The correction first identified 40 coordinate-level stereochemical errors in
three pairs. Graph-guided local reflection corrected all 40 candidates to the
query graph; every corrected candidate passed the 19 physical checks and the
repair worker never read a native pose. Reselection then recalculated chemical
graph distances for all 220 groups and evaluated 660 rule-specific selected
endpoints from the same 1,100 corrected candidates.

Checks cover atom renumbering, aromatic and resonance symmetry, nonzero
receptor-frame translation, bounded timeout, failure cleanup and hash-sensitive
recovery. Observed candidate comparisons include up to 48 complete atom maps.

`python -m routealign.reproduce` verifies the corrected seven-stratum primary
analysis and three-stratum same-candidate selector ablation. The molecular
correction is implemented in `routealign.stereo_repair`, with a synthetic
coordinate-inversion test in `tests/test_correspondence.py`. Statistical
intervals use 20,000 hierarchical resamples; family tests and Holm adjustments
are recomputed. Independent external-model integrity review is unavailable in
this correction run; execution checks are not presented as such a review.

The pre-repair audit, full correction record and compact resolved status remain
separate frozen files so that the implementation error and its resolution can be
inspected independently.
