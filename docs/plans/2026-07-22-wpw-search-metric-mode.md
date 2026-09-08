# WPW Search-space Metric-mode Diagnostic Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

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

## Execution result

Implemented in `e7ca9ca`; focused generalized-algebra, matrix-free, fixed-H,
and full MPI build checks pass, and repeated review found no blocker.  Fresh
unpreconditioned B=6 run `20260722_task16_search_metric_no_precondition_b6`
reproduced the prior residual history and `info=40` boundary.

The 480-column search metric becomes increasingly rank-deficient: effective
rank is 437 at iteration 32, 292 at 96, and 213 at 160, with the smallest
retained eigenvalue held just above the `1E-10` relative cutoff.  Nevertheless,
the residual is not concentrated in discarded directions.  At iteration 160,
discarded fractions are `1.0231E-03` occupied and `1.6596E-03` extra; the
lowest retained decade contains `9.3282E-03` occupied and `3.6089E-02` extra.
Thus more than 99%/96% of occupied/extra residual weight remains outside the
discarded-plus-lowest-retained metric modes.  Search-space redundancy exists
but cannot account for the residual plateau.  The evidence points next to the
LOBPCG recurrence/restart update in retained directions, not to a deficient
physical occupied-W span.  No checkpoint, manifest, or RT state was published.
