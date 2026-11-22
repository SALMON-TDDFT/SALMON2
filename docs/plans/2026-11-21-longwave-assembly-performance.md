# Longwave assembly performance implementation plan

**Goal:** Reduce the unchanged Gamma mini-cell quadrature and repeated box-form work while preserving optical head/wings and all convergence/finite checks.

**Measured starting point:** Nr12^3, Nk=Nq4^3, box6^3, 3Ha, 8MPI1OMP: dense661.2403s versus optical-real-space hybrid195.7473s, W error2.22e-16. Hybrid response147.8594s (finite-q solve106.0541s, optical19.1401s), assembly47.7996s. Existing baseline and frozen benchmark executables remain unchanged.

1. Add a direct finite-grid form-factor regression and matched-resource microbenchmark against the exact old per-box arithmetic. Check anisotropic grids, fractional/negative/zero/Nyquist-like coordinates and boxes1,3,4,6,12. Demonstrate the old implementation fails the local performance target before changing it; keep hardware speed thresholds out of routine correctness CI.
2. Cache the three one-dimensional box sums within each screen_forms call, then assemble Cartesian products. Preserve the original summation and multiplication/division order. No persistent geometry cache, altered physics, threshold changes or approximate phases.
3. Contract screened*direction once per angular direction before the radial loop. Validate the resulting Gamma point/average against a frozen old implementation with nonzero wings and multiple boxes, in addition to analytic one-electron/local-field/IEEE tests. If numerical equality or speed does not hold, retain the narrower proven improvement.
4. Where justified by measurement, collect weighted quadrature vectors into bounded BLAS panels, keeping the same quadrature nodes/order, convergence check and final Hermitian/full-matrix validations. This is optional after the first measurements.
5. Run audited native and full SALMON builds, scalar/SOC/dense-pump/IEEE regressions, and actual Si Nk4 at the same8MPI1OMP with frozen-GS inputs. Compare full W and optical epsilon; record phase and whole-W speed. Do not collect GS payloads/full GS output. Save results in canonical calc2026.

Next candidate after this bounded change: finite-q static RHS batching to reuse target rotations; MPI box redistribution is not the bottleneck for this particular Nk4/box6 case (about0.05s).

## Implemented scope

Retained the narrow axis-sum reuse in step 2. The old implementation failed the
optional local 2x target (1.002x); the cached implementation passes with errors
below 2.5e-16 across all box counts. Step 3 reassociation and optional BLAS panels
are deferred: this change leaves the entire longwave algorithm and its operation
order intact, so the full Si comparison tests this single optimization directly.
Review requested finite checks in the new regression; these now use the existing
scalar-IEEE-compliant `screen_finite` helper. Fujitsu compilation is untested.
