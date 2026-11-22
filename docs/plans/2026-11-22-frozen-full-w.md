# Projection-free full-q snapshot W implementation

**Goal:** Implement the user-authorized full-q occupied-orbital W update and
compare fixed/updated W for the same Si weak-pulse test.

**Architecture:** A native Fortran coordinator streams occupied k sectors and
delegates real-space H actions to their owners, reusing frozen_solve for all
charge vertices within the cutoff and both eta values. It retains the entire
coherent occupied density in the commutator. This is adiabatic snapshot RPA,
not a two-time nonstationary RPA calculation. No instantaneous eigenstates or
physical conduction manifold are constructed.

Retain all local-field G/G' entries. Check and report anti-Hermitian and q/-q
reality defects before taking the reciprocal Hermitian static part, requiring
each relative defect <= hsex_w_eta_tol. Use screen_assemble with conversion
-chi_physical*Nk*volume/2 (legacy positive-response convention). Reject any
unconverged response, nonpositive dielectric, excessive eta dependence or
nonreal/nonreciprocal final W. Preserve existing W until complete validation
and atomic checkpoint publication. Reuse the finite-q minimum-shell mini-BZ
closure for the experimental insulating head; do not claim a new optical
head/wing or metallic Drude closure. Initially require Nq=Nk, powers of two,
scalar closed-shell occupied-only RT and CPU k-only MPI. Other W modes remain
unchanged. Restart of the new mode is not claimed by this test.

Use hsex_w_solver='frozen', hsex_frozen_rank/tol and existing update schedule.
Keep solver='auto' single-q diagnostics backward compatible. Save full W as
HXT1 with unchanged TB-mBJ generating metadata, without fabricated eigenstate
occupations. Refresh endpoint HSEX immediately after update and account for
the exchange-energy jump. No production Python orchestration.

1. Add a failing small MPI full-W test (all q, local fields, uniform charge
   gauge mode, W reciprocity, rank cap and remote action failures).
2. Implement full response/assembly coordinator and run audited native tests.
3. Add opt-in adapter/update/publication/input checks and official manual.
4. Audit full build and all application objects; verify old off-path regression.
5. Generate a common equilibrium initial W using this new closure; compare
   identical 4-occupied-orbital fixed/1fs-update trajectories for Si Nr8^3,
   Nk=Nq2^3, box2^3, cutoff1Ha, eta0.01Ha, 4fs weak pulse, 2fs propagation,
   8 MPI x 1 OMP. Verify actual publication, refresh, orbitals, and pre-update
   equality. Measure full-W time, residual and symmetry/eta errors.
6. Save inputs, derived diagnostics, plots and audits to canonical calc2026;
   no GS data collection and no existing remote/fixed binary modifications.

## Exact residual bound used by the full-q path

The first local pilot showed residual-factor QR dominating each many-vertex
solve, and was stopped with its executable and logs preserved. The full-q
path now selects a conservative spectral-frame residual bound. Let L,R be
the Krylov frames rotated by the small projected eigensystems, LS=L^H U,
RS=V^H R, Uout=U-L LS, Vout=V-R RS^H. The exact source defect is
Uout S V^H + L LS S Vout^H. Include both H-frame defects HL-L eL and HR-R eR
and the spectral denominator defect, bound each product by Frobenius/operator
norm inequalities, and retain the previous roundoff allowance. Frame operator
norms are bounded by sqrt(1+||L^H L-I||F) and its right analogue. This avoids
large source-dependent basis rotations and QR; it neither truncates the
source nor assumes stationarity. The default single-q exact norm remains.
Independent full/capped-rank residual tests and the analytic full-W model
verify this upper bound; capacity exhaustion remains a rejection.

## Completed validation

The common initial-W preparation, both1600-step RT trajectories and both
zero-field controls completed. All5 full-q computations passed the unchanged
1e-6 residual-bound criterion and the static/eta/reciprocity checks. Every
update used4 occupied orbitals without instantaneous-state projection. The
new-mode-off trajectory agrees bitwise with the previous validated binary.
GNU HSE-on/off builds and MPI1/4 tests passed, including empty k/pair owners,
analytic full-W normalization, NaN/Infinity and failure rejection. Independent
review found no remaining P1/P2 issue. Fujitsu, accelerators and restart are
untested. The finite-q-shell optical closure and coordinator scaling remain
experimental. Results and initial-W convention differences are recorded in
canonical calc2026/si-frozen-full-update-20261122/report.md.
