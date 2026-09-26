# DC HSE convergence investigation

User approved convergence-first work on MPI8 Si64 8x1x1, full WF support,
no pair pruning, one simulation at a time. Preserve compiler workaround and
BLAS1/OMP2. Do not infer convergence from mixed-density changes alone.

1. Opt-in export of Psi/Hpsi before/after subspace diagonalization, CG and
   Gram-Schmidt, using a fixed local potential and ACE throughout the solve.
   Verify absence on current executable (RED), then add guarded diagnostics.
2. Run a short unchanged-physics diagnostic. Measure Rayleigh residual,
   orthogonality and occupied-subspace residual independently. Identify the
   first stage that loses accuracy before changing a solver or SCF controls.
3. Repair a demonstrated solver defect if found, with a regression. If the
   fixed-H solver is sound, test nested density SCF with frozen ACE, separate
   tolerances and fully synchronized fragment occupation/chemical potential.
   PC-DIIS is a later alternative, not a simultaneous change.
4. Run relevant complex LCFO/HSE/LAPACK regressions serially. Repeat the
   full-support MPI8 baseline, measure true density and eigen residuals;
   require convergence before any locality cutoff study. Record failures.

## Evidence and implemented correction

At iteration 20 in the unchanged-physics diagnostic, fragment1 occupied RMS
residuals (Ha) were .485 after subspace diagonalization, .317 after CG20,
1.636 after Gram-Schmidt and 1.386 after the SCF update. The Gram matrix after
CG had minimum eigenvalue 7.90e-5 and maximum overlap error .764. The returned
projected Hamiltonian had maximum off-diagonal magnitude 1.657 Ha. Thus small
individual CG residuals did not establish a solved orthonormal eigenbasis.

DC-HSE now performs a final fixed-ACE Ritz diagonalization after Gram-Schmidt.
Before occupation redistribution, calc_eigen_energy is evaluated with saved/
restored hse_freeze, so its usual exchange-refresh side effect is suppressed.
The global DC chemical potential sees the current Ritz energies and core norms.
Neither source-factor locality nor exchange/ACE formulas were changed.

Regression check_ritz.py fails on saved pre-fix data and passes on new data.
Independent check_occupations.py reconstructs the chemical potential using the
current Ritz spectrum and core weights: old data max occupation error2.0,
new data4.20e-12; global core electron count agrees to1e-7. Stage-output smoke
check_solver.py failed with missing outputs before instrumentation, then passed.

With CG20, corrected final-current-H residual at20 was .273Ha in fragment1.
Changing only ncg to1 reduced it to .188Ha and the after-CG overlap error from
.905 to .00268. Running CG independently for many steps lets bands approach
each other; frequent orthogonalization/Ritz diagonalization mitigates this.
The production default ncg has NOT been changed. The new benchmark input uses1.
Seven complex LCFO130/HSE422/LAPACK CTests pass serially.

Read-only review found no blocker in the fixed-exchange evaluation, exchange
refresh ordering or diagnostics scratch reuse. Stage snapshots retain occupations
BEFORE redistribution; weighted stage residuals must be labelled accordingly.
The final snapshot uses a refreshed Hamiltonian, whose eigenvalues need not
match the earlier occupation assignment before overall self-consistency.

Opt-in SALMON_HSE_SOLVER_DIAGNOSTIC=1 writes four version1 Psi/Hpsi files,
overwriting them at each solve (not suitable for long production runs).
The after_orthogonalization tag includes the final Ritz diagonalization for
DC-HSE. Other cases retain their prior solver behavior.

## Converged baseline

One MPI8/OMP2/BLAS1 run with ncg=1, ncg_init=4, simple mixing .01 reached
rho_dne<1e-7 at iteration1148 in595.927s (native convergence, no timeout).
All eight snapshots report SCF convergence. Refreshed-H occupied RMS residuals
range6.427e-9--1.100e-8Ha; maximum all-state residual2.299e-6Ha (includes empty
states); maximum orthogonality error1.655e-14. Core electron sum is
256.0000000005075, with each core32 within4e-9. Mixed density residual9.982e-8
corresponds to an unmixed simple-mixing residual9.982e-6, so the eigenstate
check is essential alongside the density stopping criterion.

No nested SCF or PC-DIIS was needed to reach this baseline. There is no
measured speedup comparison from this run. The reference input now records
these settings. Full real-space support and all pairs were retained; MLWF
optimization was deliberately limited to one iteration and is NOT converged.
The SCF state can now be used for a separate localization study, but no WF
cutoff is certified. Provenance and per-fragment scalar diagnostics are in
samples/dc_hse/si64-chain/converged-status.json. Full outputs and binary
snapshots remain in workspace work/si64-chain/ritz-cg1-reference.
