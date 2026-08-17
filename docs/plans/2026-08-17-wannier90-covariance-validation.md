# Wannier90 Covariance Validation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the premature individual-centre gate with a measured affine-generator covariance gate for the converged Wannier90 transform.

**Architecture:** Add a focused MPI-safe covariance validator that consumes one band/Wannier representation pair at a time and returns a maximum defect plus checked receipt. Production reconstructs each generator action using existing row-owned assembly, validates the W90 transform, and retains the later post-character individual-centre gate.

**Tech Stack:** Fortran 2008, MPI, BLAS, Python route contracts, existing SALMON DG MPI fixtures.

---

### Task 1: Add covariance validator RED

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Add a two-channel dense representation fixture with a covariant transform.
2. Add a noncommuting transform fixture that must reject.
3. Run `python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py` and verify the missing API fails compilation.
4. Implement the validator with finite/shape checks, rank agreement, checked allocation, `U^H D_band U-D_wann`, and maximum-defect rejection.
5. Re-run the MPI 1/2/4/8 fixture and expect PASS.

### Task 2: Replace the production gate

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Add a route RED requiring generator covariance after W90 and forbidding the immediate centre-orbit gate.
2. Preserve the spectral AMN until post-W90 validation.
3. Stream each affine generator representation, form its spectral Wannier representation, and call the covariance validator.
4. Accumulate the maximum defect/workspace and provenance; then release AMN.
5. Keep the later post-character periodic-position centre-orbit gate unchanged.
6. Run the route checker and focused MPI fixtures.

### Task 3: Verify Si64

**Files:**
- No source changes expected.

1. Build `/tmp/salmon-wpw-full`.
2. Run the Si64 case with 8 MPI ranks and one BLAS/OpenMP thread per rank.
3. Verify Wannier90 converges, generator covariance passes, and execution advances beyond the former immediate centre-orbit stop.
4. Record the next authoritative outcome and resource usage.
5. Run `git diff --check` and preserve unrelated user-owned dirty files.
