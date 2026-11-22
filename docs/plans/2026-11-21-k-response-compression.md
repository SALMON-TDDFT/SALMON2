# Finite-k shared response-space implementation plan

**Goal:** Measure mathematical compression of repeated Sternheimer solves at a fixed target k before changing production W.

**Architecture:** Opt-in equilibrium diagnostic gathers occupied source states, chooses short nonzero q points, and forms the same folded charge vertices as the production driver. A shared projected Galerkin space treats each source energy as an exact shift, without energy interpolation. Rank-revealing right-hand-side selection and residual enrichment determine its size; compare independently converged scalar CG solves and true residuals. Keep W publication disabled for this diagnostic and all default paths unchanged.

**Tech Stack:** Native Fortran, existing real-space Hpsi, MPI state collection and collective rejection, BLAS/LAPACK. Python is offline test/report tooling only.

1. Add failing analytic tests for a shifted multi-RHS response API: dependent complex RHS, noncommuting occupied projector, differing shifts, full-rank agreement, exhausted capacity and NaN/Inf rejection. Preserve the red build log.
2. Implement bounded shared response-space solver in `src/xc/hsex_response_space.f90`: double-orthogonalized pivoted residual expansion, reusable H actions, small Hermitian eigensolve and exact shifted denominators, explicit final Hpsi residual check. Count all H actions and retain failure status.
3. Add `hsex_reduced_target` (0 disables), `hsex_reduced_qcount`, `hsex_reduced_maxrank`, `hsex_reduced_tol` to namelist/default/broadcast/log/checks. Restrict to CPU full-k scalar closed-shell insulating equilibrium GS, no spatial/orbital decomposition. Diagnostic is isolated from production W publication.
4. Add fixed-k comparison with shortest nonzero q selection and production G/cell/folding conventions. Report RHS singular spectrum, retained subspace, per-q target response errors, solution/true residual errors, solve/build/validation wall time and all Hpsi calls. Include independent source energies and do not add q and -q a second time. This diagnostic is one target-k contribution, not a complete W or optical spectrum.
5. Source-audit before builds; object-audit all linked application objects. Run analytic tests, 4MPI12OMP collective/layout test and existing Si Nr12 Nk4 scratch GS target-k comparisons. Reuse GS in place; do not collect GS data or output logs. Never change remote production jobs or fixed binaries.
6. Document controls/limits and save filtered diagnostics, inputs, derived measurements, source/build audit records under canonical calc2026. Benchmark includes compression setup and independent validation separately. Report negative speedup or rank exhaustion honestly; do not enable the method by default.
