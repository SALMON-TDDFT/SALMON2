# Periodic density grid index repair

Human explicitly requested repair, validation and commit/push after the prior diagnosis. Apply the already described integer-rounding then periodic-modulo design; do not change grid-commensurability tolerance, symmetry group or physical settings. Do not resume stopped calculations.

1. Add checked production-module regressions demonstrating negative roundoff, N boundary rounding, multi-period translations, valid half-cell mapping and invalid grid rejection. Observe failures before implementation.
2. Use one private index function from both init_sym_rho and sym_rho. Check finite coordinates and preserve1e-3 commensurability criterion; wrap rounded integer with modulo into0..N-1. No shared state/new input variables.
3. Run checked array tests with OMP1/2, serial and MPI builds/runs. Audit source before builds and every linked application object before executing SALMON. Exercise real W/Se skew-cell input with original route preserved and only this focused file changed. This is bounded verification, not a restart of physical5fs campaign.
4. Independent review, manifests/provenance, commit/push same approved repository/branch. Preserve historical source/output, failed binaries and all old locks. Report untested optional paths and retain all prior results as unaccepted.

## Verification outcome

- Regression452 demonstrated the original failure before the fix. The repaired production module passes16 valid boundary/density/charge subcases with OMP1/2 and3 invalid-input checks (noncommensurate, NaN, just outside unchanged1e-3 tolerance).
- GNU serial Debug with array bounds checks and MPI Release builds: CTest covariance +450/451/452 all10/10 passed for each configuration. Source and all application-object prohibition guards passed.
- Independent original-route MPI2 build retains -fcheck=all with Release optimization. All9 actual W-z/Se-z/Se-x inputs at target0.32/0.30/0.28Angstrom completed one SCF iteration and MPI finalization with exit0 and finite energy, valid normalized k weights (48/84/114 representatives). No Fortran runtime bounds error. Existing conjugate-gradient array-temporary warnings are preserved in stderr; these were not called zero stderr.
- Initial unoptimized Debug real-input check reached its90s diagnostic timeout; evidence retained, no children remained. Optimized checked confirmation used new directories; this was bounded validation, not restarting5fs physical campaigns.
- These tests establish periodic mapping and bounded application execution, not converged GS,5fs current convergence, production performance, or recovery of old invalid outputs. No stopped physical calculation was resumed.
- Independent review found no important issue. GPU/OpenACC, ScaLAPACK, EigenExa, FFTW, LibXC and alternate #ifdef test implementation were not exercised. Source changes do not alter those feature guards or communication layout.
