# Direct Proca coefficients implementation plan

User-approved design: specify paper coefficients a2,a0; map to SALMON's existing integrator.
Continue on the requested TDCDFT branch. Keep the 4³ k grid fixed; do not perform k convergence.

Goal: accept tdcdft_a2 and tdcdft_a0 in &functional for tdcdft='proca', retaining old inputs.
Architecture: SALMON a=Axc/c has the opposite vector-potential sign to Eq. (1) of
Dewhurst et al. (PRB 111 L060302). Its equation is a2*a''+a0*a=-4*pi*j_number.
Normalize alpha=-4*pi/a2, restoring=a0/a2. a2 is dimensionless; a0 uses inverse time squared
in the selected unit system. Permit both signs of a2; require nonnegative a0/a2. a0=0 is allowed
as the massless reference. Optional existing damping remains a normalized inverse-time coefficient.
Reject simultaneous nonzero alpha/restoring and direct coefficients. Log raw and effective values.
The mapping is one-to-one for finite nonzero a2, so the existing checkpoint's normalized coefficients
already enforce coefficient compatibility; preserve its format and legacy Proca restart support.

1. Tests first: extend testsuites/142_Si_tdcdft/test_update.f90 to check the sign and oscillator
   response for both-sign Proca coefficients; observe failure before implementing the helper.
2. Add proca_coefficients to src/rt/tdcdft_lrc.f90. Add input globals, namelist, defaults,
   broadcasts, a0 unit conversion, normalization before logging, and compatibility checks in
   src/io/salmon_global.f90 and src/io/inputoutput.f90. Keep the timestep update unchanged.
3. Extend test_regression.py: equivalent normalized/direct trajectories, massless LRC equivalence,
   positive a2 support, conflicting/zero/opposite-sign/nonfinite inputs, unit conversion,
   unchanged/changed-coefficient restart; run serial and MPI regression.
4. Obtain Si parameters from the paper's supplemental table, recording exact source. If inaccessible,
   explicitly label the alpha-matched pair a2=-20*pi,a0=-0.2 as a comparison rather than paper values.
   Update sample inputs, README and SALMON-DOCS patch. Never infer exact paper a2 from a plot.
5. With verified paper parameters, run a 4³ Si impulse calculation against existing ALDA/LRC data,
   save exact input, results and caveats. No k-point convergence scan.
6. Review, appropriate builds/tests, commit on TDCDFT; no push requested.

## Reference change requested by user

Williams–Ullrich arXiv:2501.13290v2 is now the primary reference. Retain the tested
optional a2/a0 conversion, but use alpha/beta/gamma in the Si templates and documentation.
Keep alpha=0.2 and beta=0; gamma=1e-4 is an explicitly unvalidated sensitivity trial,
not a fitted Si parameter. Compare against gamma=0 and 4e-4 before interpreting pump changes.
The old Dewhurst-parameter full run was stopped at about step 10780/12000; partial raw
data remain in calculations/si_tdcdft_k4/proca and are not a completed comparison.
No instantaneous counter-force cancellation or automatic gamma prescription is added.
The 4³ k grid and existing completed ALDA/LRC results remain the comparison baseline.
