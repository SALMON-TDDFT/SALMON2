# WPW Metric-aware Correction Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Compare the current unpreconditioned LOBPCG residual correction with a fragment-local `S_block^{-1}r` correction.

**Architecture:** Add a default-off, diagonal-mutually-exclusive input control and route the existing validated metric block-Jacobi inverse through a collective production adapter for fixed-H and continuation only. Preserve the normal production route's legacy no-callback behavior and every solver/publication decision outside correction-vector construction.

**Tech Stack:** Fortran, MPI, Python source contracts, CMake.

---

### Task 1: Add the metric correction comparison route

**Files:**
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`
- Modify: `tests/dg/test_dg_wpw_production_operator_mpi.f90`
- Modify: `src/common/dg_wpw_bounded_operator.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`

1. Add RED source assertions for declaration, namelist membership, default `n`,
   broadcast, variables log, y/n validation, and rejection when both diagonal
   and metric preconditioners are `y`.
2. Add RED assertions that fixed-H and continuation select exactly one of
   diagonal, metric block, or no callback, while the normal production algebra
   route remains callback-free under default controls. Preserve retained-search
   propagation in every branch.
3. Add a RED production-adapter contract requiring collective metadata validation
   on the bounded-operator communicator before delegation on every rank.
4. Extend the production-operator two-rank MPI fixture to invoke the real adapter
   on its non-identity SPD block, verify the exact local block solve, verify
   continued validity after `set_dg_wpw_interface_lambda`, and verify rejection
   after a genuine operator epoch/layout fingerprint change.
5. Run the matrix-free and production-operator tests; expect failure because
   the control and route do not exist.
6. Add `yn_dg_wpw_metric_preconditioner='n'` to global/input plumbing and enforce
   mutual exclusion with `yn_dg_wpw_preconditioner` during argument validation.
7. Add a public common adapter with the eigensolver callback signature.  Reduce
   eigenvalue/RHS validity collectively on `op%w_schedule%comm`, require finite
   eigenvalues and matching RHS count, then make every rank delegate to
   `apply_dg_wpw_fragment_block_preconditioner`; eigenvalue values do not enter
   the block inverse.
8. Add the LCFO wrapper that calls this adapter.  Branch fixed-H and continuation
   through metric / diagonal / none.  Leave `wpw_algebra_step` callback-free and
   retain its search-history keyword.
9. Log the effective correction mode and test all configurations: default
   diagonal=y/metric=n, none n/n, metric n/y, and invalid y/y.
10. Run the matrix-free and production-operator tests; expect PASS.
11. Run `git diff --check` and commit only intended implementation/test files.

### Task 2: Verify, review, and compare B=6

**Files:**
- Modify: `docs/plans/2026-07-23-wpw-metric-correction.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run matrix-free, fixed-H, generalized-algebra, occupied-W MPI tests, full
   MPI/EigenExa build, and `git diff --check`.
2. Request code review for mutual exclusion, callback lifecycle/fingerprint,
   MPI collective ordering, route-specific legacy behavior, and default
   compatibility.
3. Resolve every Critical/Important finding and rerun verification.
4. Commit only intended source/test changes.
5. Create a fresh Task-16-equivalent B=6 run with
   `yn_dg_wpw_preconditioner='n'`,
   `yn_dg_wpw_metric_preconditioner='y'`, and retained search history.
6. Run 8 MPI ranks through the unchanged fixed-H boundary.
7. Record the fragment block dimension/condition and compare inner 32/96/160
   occupied/extra residuals, correction ratios, effective rank, Ritz boundary
   defects, final `info`, and publication state with Task 16. Interpret failure
   only as failure of the local block-Jacobi metric approximation.
8. Record the evidence and next hypothesis in this plan and the session handoff.
9. Run `git diff --check` and commit only the two result documents.
