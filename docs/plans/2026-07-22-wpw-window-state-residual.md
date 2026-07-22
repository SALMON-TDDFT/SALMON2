# WPW Window State-residual Diagnostic Implementation Plan

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
   `[DG-WPW-WINDOW-STATE-RESIDUAL]` output after raw convergence evaluation.
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
