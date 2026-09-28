# Fixed EXX radius priority correction

User-specified contract: positive exx_mlwf_radius determines the spherical
support. exx_mlwf_norm_fraction is then a diagnostic target only; if omitted
or zero the warning target is .999. Never enlarge R or reject a calculation
solely because retained norm is below the target. Without positive R, keep
adaptive quantile selection. Fraction means squared norm, not geometric volume.

Implement in both existing static SCF source paths: full-k and spatial/orbital
MPI. Preserve the existing fixed-radius static-only restriction, ambiguous-center
protection, localization convergence guards and DC canonical full-support default.
Do not claim new force/MD or fixed-radius RT support. Fixed-radius pair screening
remains unsupported (previous simultaneous controls were rejected).

1. RED: source probe checks fixed R independent of fractions .5/.999/1, loss
   reporting and periodic/MPI geometry. Input/native tests expect combined
   controls accepted and equal energy to radius-only; insufficient norm warns.
2. GREEN: optional fixed radius in distributed mask; radius-first input/routing;
   consistent SCF readiness and warnings in both source backends.
3. Verify ScaLAPACK build, MPI1/2/4 source probe, full-k and spatial SCF tests,
   existing input/pre-SCF regression, then update docs and commit.

## Verification completed

- ScaLAPACK-enabled build completed successfully.
- Radius-priority, input, PBE pre-SCF and pair-screen regression suites: 33 tests
  passed. HSE and PBEh spatial SCF both retain the requested R and warn when the
  retained squared norm falls below the diagnostic target.
- Distributed source probe passed with MPI 1/2/4. Changing the diagnostic
  fraction at fixed R produces bit-identical masked sources; SCF energy
  comparisons use the existing 2e-6 eV convergence tolerance.
- Independent review identified a full-cell-radius edge case: no truncation
  should require no converged localization gauge. Reproduced the failure, fixed
  the spatial support activation guard, and verified R=100 bohr with a deliberately
  unconverged gauge (one localization iteration).
- Fixed-radius RT/MD and fixed-radius pair screening remain unsupported.
  DC canonical sources remain full-support unless localization is explicitly
  opted into with yn_exx_dc_mlwf='y'.

Local test logs: work/radius-priority-final.log, work/radius-priority-probe.log,
work/radius-priority-build.log (relative to the parent workspace).
