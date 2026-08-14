# Overlap Tile Loop Optimization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove repeated tile allocation and transpose-copy overhead from distributed symmetry overlaps while preserving row-owned memory and numerical results.

**Architecture:** Reuse maximum-width tile buffers across the full symmetry/tile loop and form the Reduce-scatter send layout directly with BLAS MATMUL.  Keep the existing MPI ownership and reduction semantics unchanged.

**Tech Stack:** Fortran 2008, MPI, BLAS-backed MATMUL, Python route/source regression tests, CMake.

---

### Task 1: Add the structural performance regression

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1:** Add checks requiring tile buffers to be allocated outside the symmetry/tile loop and rejecting `local_tile` in the overlap assembler.

**Step 2:** Run `python3 tests/dg/check_dg_overlapping_wannier_route.py` and verify RED against the current implementation.

### Task 2: Reuse overlap tile storage

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

**Step 1:** Change the tile width to 64 and precompute the maximum local owner extent and full-width receive counts.

**Step 2:** Allocate `image_tile`, `packed_tile`, and `reduced_tile` once before the operation loop, using active slices for each tile.

**Step 3:** Replace `local_tile` plus transpose with the equivalent direct tile-major MATMUL.

**Step 4:** Update workspace accounting for the persistent reusable buffers and ensure all failure exits clean them up.

**Step 5:** Run the route checker and verify GREEN.

### Task 3: Verify numerical and MPI behavior

**Files:**
- Test: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

**Step 1:** Run the construction fixture on MPI 1/2/4/8.

**Step 2:** Build the Release target with `cmake --build /Users/otobetoshihito/SALMON-dev/verification/20260814-core-ownership-fix-build -j 1`.

**Step 3:** Run `git diff --check` and inspect the scoped diff.

**Step 4:** Commit only the overlap implementation, regression, and plan documents.

