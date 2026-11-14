# ROSHAMBO reproduction protocol

Use the author implementation at commit `916fe97f81a3fdb48e7bb85384eb597836d2b01e` from
https://github.com/molecularinformatics/roshambo.

Use its separate RDKit 2023.09.6 environment. Run seeds 101, 127, 151, 179 and 199,
requesting ten ETKDGv2 conformers per seed with MMFF94s optimization, Gaussian volume,
carbon radii, color scoring and ComboTanimoto ranking. Restore the author program's
output to the supplied template coordinate frame before endpoint evaluation.

Generation requires an explicitly selected GPU and an exclusive per-device lock.
The wrapper restricts upstream CPU pools to the task affinity mask. Endpoint
evaluation uses the common RDKit 2026.09.1 / PoseBusters 0.5.1 environment and hides
all GPUs. Native query coordinates are used only by the endpoint evaluator.
