# WPW Window State-residual Diagnostic Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Identify whether occupied or extra retained states control the fixed-H window-solver residual plateau without changing solver behavior.

**Architecture:** A pure helper summarizes per-state raw and preconditioned residual norms into occupied/extra maxima, worst indices, and response ratios. `solve_window` computes the extra diagnostic Gram only at selected iterations and prints one rank-zero record. Existing scalar convergence and failure paths remain authoritative.

**Tech Stack:** Fortran, MPI, LAPACK, Python source-contract tests, CMake.

---

### Task 1: Add the pure state-block summary

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Add RED tests for occupied/extra worst-state selection, exact ties,
   one-state extra block, invalid `nocc`, negative/nonfinite norm squares, and
   zero raw residual in the preconditioner ratio.
2. Run the focused matrix-free MPI test and verify the missing helper failure.
3. Implement the pure helper with deterministic first-index tie breaking.
4. Rerun the focused test and verify PASS.

### Task 2: Instrument solve_window without changing convergence

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`

1. Add RED source-contract assertions requiring `nocc` in `solve_window`, the
   selected-iteration predicate, a diagnostic-only preconditioned Gram, and
   `[DG-WPW-WINDOW-STATE-RESIDUAL]` output after raw convergence evaluation at
   selected iterations, including capped failure iteration 160.
2. Run both source-contract tests and verify RED.
3. Pass `nocc` through both matrix-free callers.  At selected iterations,
   compute per-state raw and preconditioned norm squares with the existing
   global Gram callback, call the helper, and print the split summary on rank
   zero.
4. Ensure the existing `residual<tol .and. orth<tol` return and `info=40`
   remain unchanged and independent of the new summaries.
5. Rerun focused tests and verify PASS.

### Task 3: Verify, review, and rerun B=6

**Files:**
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run focused matrix-free, fixed-H, and occupied-W MPI tests.
2. Run `cmake --build build-mpi-eigenexa -j4` and `git diff --check`.
3. Request code review and resolve all Critical/Important findings.
4. Commit only the diagnostic implementation/tests/plans.
5. Run a fresh eight-rank B=6 calculation.
6. Record occupied/extra residual histories, worst state IDs,
   preconditioner-response ratios, stopping point, and publication status.

## Execution result

Implemented in `920381e`; focused matrix-free, fixed-H, and occupied-W tests,
the full MPI build, and two review passes completed without a remaining
Critical or Important issue.  Fresh B=6 run
`20260722_task14_window_state_residual_b6` reproduced the structural seed
diagnostics and reached the same `info=40` fixed-H boundary.

Both blocks remain unconverged at iteration 160: occupied maximum
`4.8659E-04` at state 66 and extra maximum `1.6650E-03` at state 146.  Thus
extra states set the scalar maximum, but accepting only the occupied block
would still miss `1E-8` by about four orders of magnitude.  The diagonal
preconditioner amplifies the raw-worst residual: at iteration 160 the occupied
norm changes from `4.8659E-04` to `2.2726E-02` (ratio `46.704`), and the extra
norm from `1.6650E-03` to `5.6867E-01` (ratio `341.54`).  Ratios over the
selected history range from about 33 to 563 for occupied states and 8 to 910
for extra states.  No WPW checkpoint, manifest, or RT state was published.
