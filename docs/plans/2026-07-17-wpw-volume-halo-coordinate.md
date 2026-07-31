# WPW Piecewise-Wannier Volume Correction Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make every WPW volume and density consumer use fragment-owned Wannier rows only while preserving bounded P support, canonical face coupling, and `orthonormal_ww`.

**Architecture:** Freeze the piecewise-Wannier row-layout invariant before changing production. Remove remote W from volume accumulation, WP publication, and density consistently, while retaining support W in the face scanner and leaving the standalone halo provider tests intact.

**Tech Stack:** Fortran 2008, MPI, CMake, Python test launchers, SALMON DG/WPW production path.

---

### Task 1: Freeze the piecewise-Wannier volume invariant

**Files:**
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`
- Modify: `tests/dg/check_dg_wpw_nonowned_candidate_staging.py`
- Test: `src/gs/dc/lcfo_flux.f90`

**Steps:**

1. Require the volume accumulator W dimension to equal `size(wpw_owned_w_ids)`.
2. Require WP volume publication to enumerate owned W IDs, not support W IDs.
3. Require the density path to multiply owned-W points by owned coefficients while P remains support-distributed.
4. Require the face scanner to continue using `wpw_support_w_ids`.
5. Run both source-contract tests and confirm the new volume invariant fails before production correction.

### Task 2: Correct production volume and density row layouts

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/check_dg_wpw_production_consumer.py`
- Test: `tests/dg/check_dg_wpw_nonowned_candidate_staging.py`

**Steps:**

1. Size and populate the volume accumulator with `wpw_owned_w_ids` only.
2. Build and publish WP volume candidates as owned-W by support-P.
3. Fetch support coefficients as before, then select owned W coefficients for density and potential updates.
4. Remove production volume-W halo preparation/reachability, retaining support-W face preparation.
5. Run the two source-contract tests and require PASS.

### Task 3: Verify distributed halo and quadrature consumers

**Files:**
- Modify only if a newly exposed invariant needs coverage: `tests/dg/test_dg_wpw_volume_halo_periodic_images_mpi.f90`
- Test: `tests/dg/test_dg_wpw_rank_local_quadrature_mpi.f90`

**Steps:**

1. Run all `tests/dg/run_dg_wpw_volume_halo*.py` launchers to ensure the standalone API remains valid.
2. Run `python3 tests/dg/run_dg_wpw_rank_local_quadrature.py`.
3. Run `python3 tests/dg/run_dg_wpw_rank_local_quadrature_mpi.py`.
4. Fix only failures caused by the corrected mapping; preserve fail-closed epoch and duplicate-source behavior.

### Task 4: Recheck solver and production contracts

**Files:**
- Test: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Test: `tests/dg/check_dg_wpw_production_consumer.py`
- Test: `tests/dg/check_dg_wpw_face_trace_scanner.py`

**Steps:**

1. Run the matrix-free SCF wide-spectrum fixture and require its eigenvalue, residual, and orthogonality checks to pass.
2. Run the production-consumer and face-trace source contracts.
3. Run `git diff --check`.
4. Build with `cmake --build build-mpi-eigenexa -j4`.

### Task 5: Validate the real two-fragment smoke and clean diagnostics

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify if needed: `src/gs/dc/dg_wpw_lda_hartree.f90`
- Test input: `/tmp/wpw_gs_smoke/inputfile`

**Steps:**

1. Run `mpirun -np 2 build-mpi-eigenexa/salmon < /tmp/wpw_gs_smoke/inputfile` from `/tmp/wpw_gs_smoke`.
2. Require real-space WW charge to match the coefficient-space WW norm within quadrature tolerance.
3. Require finite density, charge, potential, energy, and algebra residual through every configured iteration.
4. Remove verbose success-path density/potential diagnostics once the invariant passes; keep concise rank-local failure-stage output.
5. Rebuild, rerun the focused tests and smoke, and finish with `git diff --check`.

No commit, push, or pull request is performed in this session.
