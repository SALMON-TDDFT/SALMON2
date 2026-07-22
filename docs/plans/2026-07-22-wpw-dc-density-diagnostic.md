# WPW DC-density Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Allow a structurally valid normalized occupied-W seed to enter fixed-H while retaining the source-to-DC density difference as a visible approximation diagnostic.

**Architecture:** Add a small pure classifier that distinguishes invalid density norms from a valid above-tolerance approximation residual.  Use it at the existing pre-fixed-H physical diagnostic; only invalid inputs remain fatal, while a valid threshold exceedance prints a warning and continues.  All source-projection and downstream fixed-H/publication gates remain unchanged.

**Tech Stack:** Fortran 2003/2008, MPI, Python source-contract tests, CMake.

---

### Task 1: Specify the density diagnostic classifier

**Files:**
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Modify: `src/gs/dc/lcfo_flux.f90`

1. Add RED source-contract assertions requiring a dedicated classifier, a
   `[DG-WPW-DC-DENSITY-WARNING]` message, and removal of the direct
   `dc_density_residual > dg_wpw_scf_residual_tolerance` fatal predicate.
2. Run `python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py` and verify the new
   assertion fails for the missing classifier.
3. Add a pure classifier returning `residual`, `warning`, and `info`.  Accept
   finite nonnegative numerator and finite positive denominator; warn only
   when the derived residual is strictly above tolerance.
4. Replace the direct fatal comparison with this classifier.  Print the
   warning on rank zero and reduce only `info` collectively.
5. Rerun the source-contract test and verify PASS.

### Task 2: Cover numerical boundaries

**Files:**
- Modify: `tests/dg/test_dg_wpw_fixed_h_relaxation.f90` if the existing fixture
  exposes the classifier; otherwise add the smallest focused Fortran fixture
  under `tests/dg/` and register it in the existing runner.

1. Add RED cases for residual below, equal to, and above tolerance; zero
   numerator; zero/negative denominator; and nonfinite inputs.
2. Run the focused test and verify failure before exposing or completing the
   classifier contract.
3. Implement only the visibility or numerical checks required by the fixture.
4. Rerun the focused test and the Python source-contract test; verify PASS.

### Task 3: Verify and run the physical gate

**Files:**
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run `python3 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`.
2. Run `cmake --build build-mpi-eigenexa -j4` and `git diff --check`.
3. Commit only the density-diagnostic implementation, tests, and these new
   plan documents; preserve unrelated worktree changes.
4. Copy the normalized B=6 input into a fresh run directory and execute with
   eight MPI ranks.
5. Confirm the log contains the DC-density warning and reaches the next real
   fixed-H boundary.  Record fixed-H residuals, charge, density diagnostics,
   exact stopping point, and whether checkpoint/RT publication occurred.
6. Update the session handoff with the fresh evidence and request code review.

## Execution result

Implemented in `6622842`; focused source-contract/MPI tests and the full MPI
build pass, and the second review found no Critical or Important issues.  The
fresh B=6 run `20260722_task13_dc_density_warning_b6` passed both warning-only
diagnostics and all structural density-seed gates.  It reached fixed-H and
then stopped the WPW route at algebra iteration 1 with `info=40`: after 160
window iterations, the maximum retained-state generalized residual was
`1.6650E-03` versus `1E-8`, while metric orthogonality was `1.3563E-14`.
No WPW checkpoint, manifest, or RT state was published; the program safely
fell back to full LCFO and exited normally.
