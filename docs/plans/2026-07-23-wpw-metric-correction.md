# WPW Metric-aware Correction Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Compare the current unpreconditioned LOBPCG residual correction with a fragment-local `S_block^{-1}r` correction.

**Architecture:** Add a default-off, diagonal-mutually-exclusive input control and route the existing validated metric block inverse through the optional eigensolver preconditioner callback. Preserve every solver and publication decision outside correction-vector construction.

**Tech Stack:** Fortran, MPI, Python source contracts, CMake.

---

### Task 1: Add the metric correction comparison route

**Files:**
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`

1. Add RED source assertions for declaration, namelist membership, default `n`,
   broadcast, variables log, y/n validation, and rejection when both diagonal
   and metric preconditioners are `y`.
2. Add RED assertions that fixed-H, continuation, and normal production algebra
   select exactly one of diagonal, metric block, or no callback, while preserving
   retained-search propagation.
3. Extend the two-rank MPI fixture with a non-identity SPD metric and an exact
   metric-inverse callback; assert the callback is exercised and returns finite
   deterministic correction behavior.
4. Run `python3 tests/dg/test_dg_wpw_matrix_free_scf.py`; expect failure because
   the control and route do not exist.
5. Add `yn_dg_wpw_metric_preconditioner='n'` to global/input plumbing and enforce
   mutual exclusion with `yn_dg_wpw_preconditioner` during argument validation.
6. Add `wpw_metric_precondition` with the existing eigensolver callback
   signature.  Validate eigenvalue/RHS dimensions, then delegate to
   `apply_dg_wpw_fragment_block_preconditioner`; do not use eigenvalues in the
   metric inverse.
7. Route all three production algebra call sites through diagonal / metric / none
   branches and retain the existing search-history keyword in every branch.
8. Run the matrix-free test; expect PASS.
9. Run `git diff --check` and commit only the five implementation/test files.

### Task 2: Verify, review, and compare B=6

**Files:**
- Modify: `docs/plans/2026-07-23-wpw-metric-correction.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run matrix-free, fixed-H, generalized-algebra, occupied-W MPI tests, full
   MPI/EigenExa build, and `git diff --check`.
2. Request code review for mutual exclusion, callback lifecycle/fingerprint,
   MPI collective ordering, all production call sites, and default compatibility.
3. Resolve every Critical/Important finding and rerun verification.
4. Commit only intended source/test changes.
5. Create a fresh Task-16-equivalent B=6 run with
   `yn_dg_wpw_preconditioner='n'`,
   `yn_dg_wpw_metric_preconditioner='y'`, and retained search history.
6. Run 8 MPI ranks through the unchanged fixed-H boundary.
7. Compare inner 32/96/160 occupied/extra residuals, correction ratios, effective
   rank, Ritz boundary defects, final `info`, and publication state with Task 16.
8. Record the evidence and next hypothesis in this plan and the session handoff.
9. Run `git diff --check` and commit only the two result documents.
