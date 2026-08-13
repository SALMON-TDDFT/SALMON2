# Overlapping-Wannier Direct Redistribution Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove post-Wannier core copies and the intermediate orbital-owned full-grid array.

**Architecture:** Keep final buffer storage canonical, expose core points through their position map, and fuse orbital transpose with center-owner redistribution using bounded MPI batches. Retain old primitives only as a numerical reference in tests.

**Tech Stack:** Fortran 2008, MPI, Python route/memory checks, CMake.

---

### Task 1: Lock production memory lifetimes

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Add failing assertions prohibiting post-buffer allocation/copy of `ow_core_values`, `ow_core_gradients`, and production calls to the two-stage transpose/redistribute pair.
2. Require a fused direct redistribution call using buffer values and core positions.
3. Run both checkers and observe the expected failure.

### Task 2: Add fused redistribution primitive

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

1. Add a focused MPI fixture comparing the fused result with the established two-stage result.
2. Run MPI 1/2/4/8 and observe the missing-routine failure.
3. Implement checked batched direct routing to center owners.
4. Run MPI 1/2/4/8 and require identical orbital IDs and values.
5. Commit the primitive and fixture.

### Task 3: Route production through buffer views

**Files:**
- Modify: `src/gs/main_dft.f90`
- Test: `tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py`
- Test: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Replace final core copy allocation with buffer/core-position inputs.
2. Call the fused redistribution primitive.
3. Adapt core-only overlap/center consumers to bounded core views or direct indexed inputs.
4. Run the static checkers and compile `main_dft.f90`.
5. Commit the production route.

### Task 4: Full verification

1. Run focused construction MPI tests on 1/2/4/8 ranks.
2. Run W90 MPI tests on 1/2/4/8 ranks.
3. Run route, memory-lifetime, obsolete-route, and `git diff --check` checks.
4. Build a clean MPI/EigenExa/Wannier90/SPGLIB binary with low build parallelism.
