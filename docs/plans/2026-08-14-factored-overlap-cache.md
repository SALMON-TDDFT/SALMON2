# Factored Point-Cogroup Overlap Cache Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove repeated point-representation assembly and per-owner reductions from the post-Wannier factored proof.

**Architecture:** Cache all row-owned point representations once in the factored validator. Refactor overlap assembly to form a full local tile and distribute reduced row blocks with one `MPI_Reduce_scatterv` per tile.

**Tech Stack:** Fortran 2008, MPI, BLAS-backed `matmul`, SALMON MPI fixtures.

---

### Task 1: Add performance-contract REDs

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Add an optional prepared-operation-count receipt to the factored validator fixture.
2. Require the count to equal the point-cogroup order for the nontrivial fixture.
3. Run `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py` and confirm compilation/assertion failure.

### Task 2: Reduce overlap collectives

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

1. Precompute row receive counts/displacements with checked default-integer products.
2. Replace owner-loop partial arrays and `MPI_Reduce` calls with a full local tile and `MPI_Reduce_scatterv`.
3. Include the full tile and receive buffers in the checked workspace receipt.
4. Run the construction fixture on MPI 1/2/4/8.

### Task 3: Cache point representations

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

1. Assemble all point operations once in a single row-owned tensor.
2. Measure identity and unitarity directly from the cached tensor.
3. Use cached left/right slices in generator relation products; continue direct assembly for expected affine maps.
4. Return the prepared-operation-count receipt and update checked peak memory accounting.
5. Run the construction fixture on MPI 1/2/4/8.

### Task 4: Full verification

1. Run `python3 tests/dg/check_dg_overlapping_wannier_route.py`.
2. Run `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py`.
3. Build `/Users/otobetoshihito/SALMON-dev/verification/20260814-core-ownership-fix-build` with `-j 1`.
4. Run `git diff --check` and commit only the implementation, focused tests, and plans.
