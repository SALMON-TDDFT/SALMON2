# Spatial Gamma Jacobi reuse implementation plan

**Goal:** Unblock large-cell MLWF localization by using the existing Gamma Jacobi optimizer in the Gamma-only spatial EXX path.
**Architecture:** Both spatial-only and spatial/orbital refresh already construct six Gamma links with matching +/- weights. Consume those links via gauge_minimize_gamma_inplace, retaining accepted-gauge transport/fallback and canonical DC bypass. Full-k legacy paths remain unchanged.
**Tech stack:** Fortran, MPI, existing gauge optimizer, ScaLAPACK build and Python regression harnesses.

User explicitly approved sharing the Gamma Jacobi path before remeasurement on 2026-09-28.

1. Add a deterministic mixed periodic Gaussian subspace probe requiring converged localized sources in pure spatial and combined orbital layouts. Demonstrate failure of the generic minimizer; additionally retain the observed 64-H2 failure at gradient 1.8626923e-6, iteration406, maxiter1000/tol1e-6.
2. Replace the two Gamma-only refresh calls with existing in-place Jacobi calls. Verify source/gauge unitarity and layout parity; retain existing exchange-action/ACE and canonical DC regressions.
3. Rebuild after the independent GS preparation batch ends, preserving its executable identity. Verify stopped 64-H2 fixture; rerun benchmark with one new binary. Preserve old partial measurements and never combine their timings with the new run.

Progress: approved design recorded. GS preparations continue without overlapping numerical test jobs or changing their executable.
- RED: mixed Gaussian probe with old objects failed at 100 iterations, gradient .012419; 64-H2 actual fixture previously failed at gradient1.8626923e-6.
- GREEN: five pure spatial/combined orbital layouts converge in six Jacobi sweeps at gradient9.85e-10 and retain orthonormality/source agreement. Exchange/ACE probe passed all 27 rank/orbital/omega configurations; 20 adaptive RT/pre-SCF/radius/pair-screen regressions passed. ScaLAPACK build succeeded.
- Read-only independent review found no blocking issues. Raw links are unused after minimization; accepted transported-gauge fallback remains safe.
- 64-H2 actual initial localization now converges in three sweeps, gradient5.1780163e-8 at the unchanged 1e-6 target. Full 16-step run in progress.
- GS preparation: seven shape payloads complete; explicitly interrupted only the final 2x2x2 preparation early so primary localization verification could proceed. Its incomplete directory remains preserved, is not imported, and will be recalculated.
