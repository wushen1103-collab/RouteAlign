# Release verification

The portable workers were checked in an isolated workspace, separate from the
saved experimental outputs. Existing environments were reused.

- Frozen statistical reproduction checked 1,389 numerical quantities across
  the hybrid and uniform-medoid analyses, each with seven strata.
- A five-candidate low-similarity case reproduced RMSD 1.1132034305068892 Å.
- A five-candidate receptor-consensus case reproduced RMSD 0.2778478559439173 Å.
- Both case values and physical-pass outcomes matched their frozen endpoints.
- Native query poses were withheld during candidate generation, repair and
  selection, then restored only for endpoint evaluation.
- The TEMPL smoke case reproduced RMSD 0.6079861979708726 Å when generated in
  its original RDKit 2024.09.6 environment and evaluated in the common endpoint
  environment. The runner rejects the mismatched generation environment.
- An injected failure returned 97 without committing a task result. An imposed
  deadline returned 124 and stopped the task. Repeating a completed task resumed
  without replacing the saved result.
- Resource tests checked a single allowed CPU, hidden GPUs, and termination of
  a child process before it could write a delayed output.

The complete original cohort is represented by frozen observations; release
validation reran representative molecular cases rather than repeating the full
experimental campaign. ROSHAMBO generation was not rerun during packaging.
Its original observations, frozen source version and controlled worker are
included in the reproduction protocol.

Machine-readable numerical checks are in [verification.json](verification.json).
