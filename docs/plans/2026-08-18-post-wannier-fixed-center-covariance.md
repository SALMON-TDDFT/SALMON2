# Post-Wannier Fixed-Center Covariance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Validate Wannier90 against the complete fixed-center DMN contract while leaving full-affine validation to the post-MLWF symmetry construction.

**Architecture:** Reuse the existing streamed covariance validator, but assemble the post-Wannier operations from the same fixed-center maps and ordering used to write DMN. Keep one dense representation live at a time and retain maximum defect/workspace receipts.

**Tech Stack:** Fortran 2008, MPI, Wannier90 checkpoint/DMN, Python route contracts.

---

### Task 1: Add a failing route contract

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Isolate the source slice from `run_dg_w90_gamma_library` through the covariance-passed diagnostic.
2. Require `fixed_center_group_order` and the matching `fixed_center_symmetry_map` column.
3. Reject `global_affine_generators` and `global_symmetry_map` in that slice.
4. Run the route checker and confirm it fails on the current global-affine loop.

### Task 2: Align the post-Wannier gate with DMN

**Files:**
- Modify: `src/gs/main_dft.f90:1719-1760`

1. Loop over `1:fixed_center_group_order`.
2. Assemble from `fixed_center_symmetry_map(:,operation:operation)`.
3. Reconstruct the target representation from the same `spectral_amn` convention used by DMN.
4. Rename diagnostics and fatal messages from generator to fixed-center covariance.
5. Run the route checker and confirm it passes.

### Task 3: Verify and rerun Si64

1. Run DMN and W90 MPI 1/2/4/8 checks.
2. Run fragment symmetry and production route checks.
3. Rebuild production and run `git diff --check`.
4. Run Si64 with MPI 8 / OMP 1 and require passage of the immediate fixed-center covariance diagnostic.
5. Continue into the downstream translation and point-cogroup gates without relaxing tolerances.

