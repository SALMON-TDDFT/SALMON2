# WPW Bounded Operator Index Cache Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Codex:** Execute with test-driven development in the existing protected worktree.

**Goal:** Remove invariant stable-ID linear searches from every bounded H/S application.

**Architecture:** Preserve global-ID sparse storage and add epoch-local integer endpoint caches built transactionally during operator initialization.  Apply kernels consume only cached positions without changing entry or summation order.

**Tech Stack:** Fortran 2008, MPI, Python source contracts, CMake.

---

### Task 1: Add the RED apply-loop contract

**Files:**
- Create: `tests/dg/check_dg_wpw_bounded_index_cache.py`

1. Extract the `apply_blocks` source body.
2. Reject any `find_id` call in that body.
3. Require six endpoint cache arrays and their initialization assignments.
4. Run and confirm failure against the existing repeated lookup loop.

### Task 2: Build and consume endpoint caches

**Files:**
- Modify: `src/common/dg_wpw_bounded_operator.f90`

1. Add six allocatable integer arrays to the bounded operator type.
2. Allocate and populate them in candidate initialization while validating endpoints.
3. Replace all six apply-loop lookups with cached array references.
4. Preserve the sparse entry traversal and accumulation order exactly.
5. Run the new contract until GREEN.

### Task 3: Verify functional equivalence

**Files:**
- Inspect: `tests/dg/test_dg_wpw_gs_bounded_apply_mpi.f90`
- Inspect: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Run the bounded H/S dense-oracle fixture.
2. Run the matrix-free wide-spectrum fixture.
3. Run production source contracts and candidate routing tests.
4. Run the full MPI/EigenExa build and `git diff --check`.

### Task 4: Re-run full-basis Si64 preflight

1. Use a new dedicated directory and the manifest-pinned atom/pseudo hashes.
2. Run 8 MPI ranks and 1 OpenMP thread with one WPW SCF iteration.
3. Require finite DC convergence, successful operator build, no rank collapse, and a bounded first algebra-step duration.
4. Review local rank diagnostics before deciding whether production GS may start.

No commit, push, or pull request is allowed.
