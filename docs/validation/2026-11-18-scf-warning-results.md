# Advisory SCF non-convergence (2026-11-18)

DFT status reports with finite energy/residual, positive tolerance, nonnegative
residual, positive iteration count and a valid completion flag now return even
when the SCF threshold is not met. RuntimeWarning includes iteration count,
residual and tolerance. `converged` is false unless flag=1 and residual<tolerance.
Completed execution is still `status=complete`; this is not a convergence claim.
The manifest retains the convergence flag and measured stopping quantities.
Legacy completed manifests without the new field remain readable; status records
provide the authoritative convergence value. New metadata must agree with them.

Malformed/missing reports, nonfinite values, missing restart data, nonzero native
exit or incomplete execution remain errors. No normalization, tolerance changes,
extra SCF or native formula changes are introduced.

Verification on this Apple Silicon Mac:

- Source suite: 581 tests passed, one optional test skipped.
- Unit regression: finite unconverged results survive save/reopen, force access,
  warm start and the TDDFT bridge. Missing/invalid reports still fail.
- Actual bundled Si starter input with max_iterations=1 and tolerance=1e-12:
  converged=false, iterations=1, residual=1.09383 (display-rounded); finite forces.
  Reopened the saved result and performed two TDDFT steps with zero external A;
  current remained finite. This verifies continuation, not physical convergence.
- Existing explicitly built PZ executable used without rebuilding or modification.
  Script and report: /private/tmp/test-unconverged-scf.py and
  /private/tmp/salmonpy-unconverged-scf-20261118/report.json.
- No GitHub Actions. New native MPI tests were not required for this Python-side
  serial DFT result policy. Pending fast-PR asset MPI admission remains separate.
