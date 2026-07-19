# WPW Candidate Routing Performance Implementation Plan

> **For Codex:** Execute task-by-task with test-driven development.

**Goal:** Replace the O(M^2) WPW candidate routing sort with deterministic O(M log M) index sorting and reusable orders.

**Architecture:** Candidate records remain immutable after MPI routing. Stable mergesort orders integer indices using the existing WP/PP comparator, and validation plus aggregation reuse those orders. MPI ownership and payload contracts do not change.

**Tech Stack:** Fortran 2008, MPI, Python test runners, CMake.

---

### Task 1: Add a performance regression

**Files:**
- Modify: `tests/dg/test_dg_wpw_candidate_halo_mpi.f90`
- Modify: `tests/dg/run_dg_wpw_candidate_halo_mpi.py`

1. Add a reverse-ordered, duplicate-free candidate set large enough to expose quadratic sorting.
2. Require the route to preserve deterministic aggregation and complete under a bounded test timeout.
3. Run the focused test and verify it fails by timeout with the insertion sort.

### Task 2: Implement stable index mergesort

**Files:**
- Modify: `src/rt/dg/rt_dg_wpw_candidate_halo.f90`

1. Add index initialization and stable bottom-up mergesort helpers.
2. Compare candidate records through indices with the existing ordering keys.
3. Change duplicate validation to consume an index order.
4. Change publication to consume precomputed WP and PP orders.
5. Remove record-copy insertion sorts.
6. Run the focused regression and verify it passes.

### Task 3: Add bounded diagnostics

**Files:**
- Modify: `src/rt/dg/rt_dg_wpw_candidate_halo.f90`
- Test: `tests/dg/check_dg_wpw_candidate_route_performance.py`

1. Write a source-contract test for candidate count and elapsed route timing.
2. Verify RED.
3. Add root-rank route count/timing diagnostics using `MPI_Wtime`.
4. Verify GREEN.

### Task 4: Regression verification

1. Run candidate halo, face scanner/binding, quadrature, matrix-free SCF, and production input tests.
2. Build `build-mpi-eigenexa` with four jobs.
3. Run `git diff --check`.
4. Run a fresh reduced Si64 preflight and confirm the route finishes.
5. Review findings before authorizing another full GS.

No commit, push, or PR is permitted for this plan.
