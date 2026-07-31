# Overlapping-Wannier Coefficient RT Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Propagate coefficients in the accepted fixed nonorthogonal overlapping-Wannier basis without entering any legacy RT route.

**Architecture:** Extend the route checkpoint to V3 with row-owned field-free and observable matrices. Use a midpoint generalized-Hermitian eigendecomposition and exponential update, and publish a separate provenance-locked RT restart.

**Tech Stack:** Fortran 2008, MPI, LAPACK focused solver, Python source contracts, CMake.

---

### Task 1: Add metric-propagation RED fixtures

**Files:**

- Create: `tests/dg/test_rt_dg_overlapping_wannier_mpi.f90`
- Create: `tests/dg/run_rt_dg_overlapping_wannier_mpi.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Steps:**

1. Define the wished-for fixed-basis propagation API.
2. Test S-norm and field-free energy conservation.
3. Test phase and retained-basis rotation covariance.
4. Test position/velocity coupling and restart-split identity.
5. Test collective provenance and stale-generation rejection.
6. Run on 1, 2, 4, and 8 ranks and observe compile RED.

### Task 2: Implement the coefficient propagator

**Files:**

- Create: `src/rt/dg/rt_dg_overlapping_wannier.f90`
- Modify: `src/rt/CMakeLists.txt`

**Steps:**

1. Validate dimensions, finiteness, Hermiticity, generations, and
   fingerprints collectively.
2. Assemble the midpoint Hamiltonian from validated observable matrices.
3. Implement the generalized-eigenvalue exponential step with LAPACK,
   caching the field-free eigensystem and re-diagonalizing a
   time-dependent midpoint Hamiltonian.
4. Add norm, energy, covariance, and restart state diagnostics.
5. Run Task 1 fixture to GREEN on 1, 2, 4, and 8 ranks.

### Task 3: Publish and load checkpoint V3

**Files:**

- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`

**Steps:**

1. Add row-owned `H0/X/V` payload and matrix fingerprints to V3.
2. Make V2 explicitly non-reusable for coefficient RT.
3. Populate V3 only after all Task 9 acceptance gates pass.
4. Test round trip, corruption, incomplete ownership, and stale
   provenance on 1, 2, 4, and 8 ranks.

### Task 4: Isolate the production RT route

**Files:**

- Modify: `src/rt/main_tddft.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Steps:**

1. Add one explicit coefficient-RT input switch.
2. Require an accepted V3 overlapping-Wannier checkpoint.
3. Call only the coefficient propagator and return.
4. Reject basis update, tail escape, stale checkpoint, and all forbidden
   fallback routes.
5. Run source contracts and focused fixtures.

### Task 5: Verification and review

**Steps:**

1. Run every overlapping-Wannier GS and RT fixture on required MPI ranks.
2. Run `git diff --check`.
3. Perform specification and code-quality reviews.
4. Resolve every Critical/Important finding.
5. Run the clean parent-prerequisite overlay build.
6. Commit only after all gates pass.
