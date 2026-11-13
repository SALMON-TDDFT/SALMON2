# ASE → GS → bands / TDDFT validation — 2026-11-13

Environment: macOS arm64, CPython 3.14, ASE 3.29.0, GNU Fortran and LP64 OpenBLAS,
with a separately built Libxc-enabled SALMON. Integration uses PZ, not Libxc.

Input: ASE bulk Si diamond, cubic 8-atom cell, lattice 10.26 bohr, 12³ real grid,
2³ k mesh, 32 electrons and 32 orbitals; repository exercise_06 scalar KY potential.
This is a coarse workflow test, not a converged physical prediction. The KY file
has no independently verified generating-XC provenance; it is used consistently
for native/wrapped comparisons.

- GS converged after 36 native Miter updates, rho_dne threshold 1e-9.
- Total energy: -31.244772732850297 Hartree.
- Independent fresh direct SALMON execution of generated input/potential:
  energy and iteration count exactly equal.
- Reopened result creates a converged 4-point band calculation.
- Two independently created TDDFT states advance for four steps, dt=0.02 a.u.,
  with opposite A/c histories. A third fresh state exactly reproduces the first's J.
- All saved GS file hashes remain unchanged after both consumers.
- One-iteration deliberately unconverged GS is rejected despite native exit success;
  failure manifest and native logs are retained.
- Result directory relocation succeeds. A corrupted occupation.bin is rejected.
- Unit tests also reject unexpected new native inputs and verify model/potential
  copies, ASE unsupported spin/PBC and independence from caller mutation.

Native change only writes the preexisting SCF convergence flag, residual, threshold,
Miter and energy to a versioned status report. No solver step or convergence criterion
was changed. Binary SHA256 is recorded in integration validation.json.

Reproduction from checkout with active interpreter and optional ASE dependencies:

```
python tools/test_ase_ground_state.py SALMON_EXECUTABLE SI_RPS_DAT NEW_DIRECTORY
```

External executable requires source/object audits and the new GS status report.
No assertion of native GS portability to untested operating systems or of arbitrary
checkpoint compatibility follows from the Python/wheel tests.

Source and repaired, separately installed Mac wheel: 61 tests passed with ASE.
Without ASE: 56 tests passed and five optional ASE tests skipped.
