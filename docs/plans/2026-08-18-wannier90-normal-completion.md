# Wannier90 Normal Completion Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Accept a normally completed Wannier90 run at `num_iter` while preserving SALMON's numerical and symmetry validation.

**Architecture:** Keep the convergence-log parser as the boundary between Wannier90 execution and SALMON result validation. Change only its success semantics: normal completion plus final state is sufficient, while the early-convergence iteration remains a diagnostic. The existing result validator continues unchanged.

**Tech Stack:** Fortran 2008, MPI fixture, Python route checker.

---

### Task 1: Specify normal-completion behavior

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Change the exhausted-iteration fixture to include Wannier90's normal completion banner and require success with `iterations == maximum_iterations`.
2. Add a truncated-log fixture without the normal completion banner and require rejection.
3. Update the route assertion to require normal-completion validation rather than early-convergence-only validation.
4. Run `python3 tests/dg/check_dg_overlapping_wannier_route.py` and the focused MPI runner; verify RED because the parser still rejects iteration exhaustion.

### Task 2: Implement the minimal parser change

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Parse the normal Wannier90 completion banner.
2. Require both `Final State` and normal completion.
3. If the early-convergence marker is absent, set the diagnostic iteration to `maximum_iterations` rather than returning failure.
4. Do not change `validate_dg_w90_result` or downstream symmetry validation.
5. Re-run the focused route and MPI tests and verify GREEN on MPI 1/2/4/8.

### Task 3: Verify integration

**Files:**
- No additional source changes expected.

1. Run `git diff --check`.
2. Rebuild the production executable incrementally.
3. Run the focused SAWF/DMN format check.
4. Inspect the diff to ensure unrelated dirty worktree files are not staged.
5. Commit only the design, plan, parser, and focused tests.
