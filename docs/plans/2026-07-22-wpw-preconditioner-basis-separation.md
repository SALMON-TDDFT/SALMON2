# WPW Preconditioner versus Basis-conditioning Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Determine whether fixed-H stagnation is driven by the diagonal preconditioner or by W/P basis conditioning.

**Architecture:** Add a default-on comparison control that conditionally passes the existing optional preconditioner callback. Add pure denominator-statistics logic and selected-iteration diagnostic state so production can correlate residual amplification with W/P diagonal denominators without changing numerical values.

**Tech Stack:** Fortran, MPI, Python source contracts, CMake.

---

### Task 1: Add the explicit comparison control

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`

1. Add RED checks for a default `'y'` input control, namelist, broadcast, log,
   validation, and conditional callback omission.
2. Verify RED, implement the minimal control, and verify GREEN.
3. Preserve the current callback path bitwise when the control is `'y'`.

### Task 2: Add denominator statistics

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add RED pure-helper tests for signs, floor hits, scale-normalized minimum,
   inverse maximum, empty block, and nonfinite input.
2. Add RED source contracts for selected-iteration diagnostic state and W/P
   `[DG-WPW-PRECONDITIONER-DENOMINATOR]` records for occupied/extra worst IDs.
3. Implement the diagnostic transaction and pure statistics without changing
   `zw=-rw/denominator` or `zp=-rp/denominator`.
4. Verify focused tests are GREEN.

### Task 3: Verify, review, and run the comparison

**Files:**
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run matrix-free, fixed-H, occupied-W MPI tests, full MPI build, and
   `git diff --check`.
2. Request review and resolve all Critical/Important findings.
3. Commit only the intended implementation, tests, and plans.
4. Run matched fresh B=6 cases with the switch `'y'` and `'n'`.
5. Compare occupied/extra histories, exact stopping boundary, denominator
   statistics, and publication status; record the conclusion in the handoff.

## Task 1 comparison result

The default-on comparison control was implemented in `a59c5f9`, passed the
focused MPI test and full build, and was reviewed without a blocker.  Fresh
unpreconditioned B=6 run `20260722_task15_no_precondition_b6` used the same
normalized seed and reached the same 160-iteration `info=40` boundary.

Removing the diagonal preconditioner materially improves convergence.  At
iteration 96, occupied/extra maxima are `3.6079E-04`/`1.1093E-03`, versus
`1.0673E-03`/`3.4053E-03` with preconditioning.  At iteration 160 they are
`2.5607E-04`/`2.0289E-04`, versus `4.8659E-04`/`1.6650E-03`.  Thus the current
diagonal preconditioner is harmful, especially for extra states.  However,
the unpreconditioned residuals remain about four orders above `1E-8` and are
nonmonotone after iteration 112, so removing it does not establish that the
W/P basis is well-conditioned or that the LOBPCG search update is adequate.
No WPW checkpoint, manifest, or RT state was published.
