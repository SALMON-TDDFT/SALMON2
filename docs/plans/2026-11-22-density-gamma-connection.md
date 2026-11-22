# Density-response Gamma-cell connection implementation plan

**Goal:** Continue the user-approved density-response W approach with an opt-in finite-q Gamma closure for the frozen occupied-only path.

**Architecture:** Existing stationary optical head/wings and the legacy shell result remain reproducible. Add `hsex_frozen_gamma='q2'` (default `'shell'`). At every nonzero q, extrapolate density chi0 from eta, eta/2, eta/4 before dielectric inversion, then verify the W change. Reuse the validated full-matrix `screen_q2_head` nested fit and spherical-cell convention. Retain only its small-q stencil on root, not a second full-Nq kernel. Replace the complete Gamma block once, without adding the old scalar head. Do not use a current response, force a pumped Drude cancellation, or introduce empty/projection states.

**Tech Stack:** Fortran, existing SALMON MPI wrappers, LAPACK, native focused tests.

1. Add native regression tests for stencil selection, insulating/metallic full-matrix cell moments, insufficient/ill-conditioned q data and nonfinite values. Add driver early rejection tests before any Hamiltonian call.
2. Implement the stencil helper and frozen-full opt-in. Apply the existing finite-q eta/positivity/Hermiticity checks before fitting. Bound the sample storage by the nested stencil and fail before publication on allocation or fit errors. Keep all existing entry points compatible.
3. Wire namelist default/broadcast/case normalization/logging/validation and SALMON-DOCS patch. Reject incompatible stationary-optical and q2 selections. Require usable q sampling at input/preflight; the existing power-of-two frozen grid therefore needs at least 8^3.
4. Run native serial/MPI tests with empty owners, old-mode regression checks, source/object audits and HSE-on/off builds. Run the small existing Si case to verify that insufficient q data are rejected before an attempted update; reuse the established stationary Si W comparison as the unchanged-path control.

Acceptance: correct Gamma replacement and no double counting in the tested discrete-q/spherical convention; eta extrapolation precedes q fitting; failed checks cannot publish W. This does not prove convergence of an actual nonstationary Si spectrum. The available small Si snapshots are Nk=2^3 and cannot validate a two-shell limit. Actual q/eta and grid convergence must be reported separately, and the option stays off by default. No remote production jobs or fixed binaries are modified. Collected records go to canonical calc2026 without GS payloads.

## Completed validation

- New `density_gamma_stencil` and transactional `density_gamma_replace` helpers retain only selected small-q matrices on root and keep final pair ownership. The input uses `hsex_w_head_tol` for q fitting and `hsex_w_eta_tol` for finite-q eta extrapolation/change gates.
- Native 1/4-MPI complete Nk=Nq8^3 zero-density-response calculations agree with an independent bare-Coulomb and analytic spherical-cell reference: Gamma contribution relative error 0.4564%, within the 5% criterion; no double scalar head is added.
- Polynomial insulating/metallic matrix moments, odd wings, empty owners, q shortage, failed fits, negative coefficients, remote NaN/Infinity and overlapping layouts pass their accept/reject tests. Existing frozen/optical failures and successes remain checked.
- An exact paired-source cancellation previously rejected by the compressed-source guard now returns a zero response without H calls. The change accepts only original identical factors with opposite coefficients; the old unresolved-near-cancellation and resonant-roundoff rejection tests remain passing.
- OpenBLAS and Accelerate native tests, HSE-enabled/disabled GNU application builds, source/all-application-object PZHEEV audits, scalar IEEE controls, and 132-column/diff checks passed. Source-only independent review reported no important findings. Fujitsu/accelerator compilation is untested.
- Real Si Nr8^3,Nk=Nq2^3,boxes2,cutoff1 Ha,8 MPI: the stationary optical W control completed in426.90 s; relative W difference2.49334e-16 and maximum absolute difference2.77556e-17 Ha from the earlier accepted control. q2 small-mesh and conflicting-input cases reject before W publication, including uppercase-input normalization.
- Actual nonstationary Si at Nk>=8 and its q/eta/grid convergence remain untested. The experimental option stays off by default. The root fit's selected-matrix workspace still grows with box count; no large-box or full large-k performance claim is made.
