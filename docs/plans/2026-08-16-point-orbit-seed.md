# Point-Orbit Seed Preservation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Construct every point-cogroup orbit vector from one unchanged normalized seed so the Si64 periodic-center gauge can form complete orbit blocks.

**Architecture:** Retain the existing row-owned and small-matrix data structures. Split one mixed point loop into a pure orbit-generation loop followed by the existing clustering and orthogonalization loop.

**Tech Stack:** Fortran 2008, MPI, LAPACK, existing Wannier90 MPI fixture.

---

### Task 1: Add a cyclic three-operation regression

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`

**Step 1:** Add a C3 point representation in a three-dimensional internal space, with three distinguishable periodic centers, and call `jointly_canonicalize_dg_sector_periodic_position_gauge`.

**Step 2:** Assert successful construction, finite defect, and three complete center entries.

**Step 3:** Run the focused Wannier90 MPI fixture and verify the test fails on the mixed seed/residual loop.

**Step 4:** Commit the RED.

### Task 2: Separate orbit generation from orbit processing

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`

**Step 1:** Before the processing loop, fill `orbit_vectors(:,1:npoint)` from the unchanged normalized `orbit_residual` seed.

**Step 2:** Run the existing center/cluster/cover logic over the stored vectors without using a mutated residual as an operation input.

**Step 3:** Run the Wannier90 MPI fixture on 1, 2, 4, and 8 ranks, both route checks, `git diff --check`, and the production build.

**Step 4:** Commit the minimal fix.

### Task 3: Re-run Si64

**Files:**
- Create: a new verification directory outside the worktree

**Step 1:** Launch the same Si64 input on 8 MPI ranks with all thread counts fixed at one.

**Step 2:** Confirm Wannier90 explicitly converges and the joint periodic-center gauge passes the former orbit-block failure.

**Step 3:** Continue through the first subsequent production checkpoint and report the new outcome, runtime, and RSS/rank.

