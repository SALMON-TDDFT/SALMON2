# WPW Search-space Metric-mode Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Measure whether fixed-H residuals are trapped in discarded or near-cutoff LOBPCG metric modes.

**Architecture:** Add a pure/testable reduced-metric spectrum helper in `dg_generalized_algebra`. At selected `solve_window` iterations, assemble `Z†R`, call the helper on the existing `reduced_s`, and print occupied/extra mode weights without feeding them into the solver.

**Tech Stack:** Fortran, MPI, LAPACK, Python/Fortran focused tests, CMake.

---

### Task 1: Add the metric-mode classifier

**Files:**
- Modify: `src/common/dg_generalized_algebra.f90`
- Modify: `tests/dg/test_dg_generalized_algebra.f90`

1. Add RED fixtures with known diagonal metric spectra and residual overlaps.
2. Verify discarded and lowest-retained-decade fractions for occupied/extra
   columns, including zero overlap, no discarded modes, and invalid input.
3. Implement the helper using the same Hermitian eigensolver and cutoff rule.
4. Verify focused algebra tests pass.

### Task 2: Instrument selected solve-window iterations

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add RED contracts for selected-only `gram(z,r)`, helper invocation, and
   `[DG-WPW-SEARCH-METRIC-MODE]` output before the reduced eigensolve.
2. Assemble the overlap only on globally uniform selected iterations.
3. Keep all existing reduced-eigensolver inputs and convergence decisions
   unchanged.
4. Run focused tests and full build.

### Task 3: Review and matched B=6 runs

1. Request code review and resolve all blockers.
2. Run fresh preconditioned and unpreconditioned B=6 cases.
3. Record metric ranks and occupied/extra discarded/near-cutoff fractions at
   iterations 1, 32, 96, and 160 in the session handoff.
