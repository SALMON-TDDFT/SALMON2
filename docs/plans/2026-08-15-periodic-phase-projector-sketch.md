# Periodic-Phase Projector Sketch Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the O(N^2) periodic-phase alignment fingerprint communication with an O(Nm) deterministic projector sketch.

**Architecture:** Compute a fixed number of `P v = A(A^H v)` probes without materializing the global projector. Hash the globally ordered probe outputs using guarded tolerance quantization, preserving target-frame and MPI-decomposition invariance.

**Tech Stack:** Fortran 2008, MPI, Python source-route tests, focused MPI fixtures.

---

### Task 1: Add the O(N^2) regression guard

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Add a focused source check for `align_dg_w90_character_sectors_by_periodic_phase` that rejects `projector_row(global_row_count)` and a global-length projector-row `MPI_Allreduce`.
2. Run `python3 tests/dg/check_dg_overlapping_wannier_route.py` and verify it fails on the existing implementation for the intended reason.

### Task 2: Implement the fixed-probe projector fingerprint

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Replace `projector_row(:)` with length-`m` sketch workspace and scalar streamed output.
2. For each deterministic probe, accumulate the local `A^H v`, reduce it once across MPI, and evaluate each local `P v` row.
3. Stream rows in global row-ID order, quantize real and imaginary parts after finite/range checks, and fold them into the fingerprint.
4. Collectively reject numerical or MPI failures and update checked workspace accounting.
5. Run the route check and verify it passes.

### Task 3: Verify numerical invariance and integration

**Files:**
- Modify only if needed: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify only if needed: `tests/dg/run_dg_overlapping_wannier_w90_mpi.py`

1. Ensure the fixture directly compares fingerprints under target-frame rotation.
2. Run `python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py` for MPI 1/2/4/8.
3. Run the route and obsolete-route checks.
4. Run `git diff --check`.
5. Build Release with one build job.

