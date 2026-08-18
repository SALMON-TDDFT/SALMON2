# DG Wannier Iteration Setting Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the DG overlapping-Wannier route honor SALMON's existing `wannier_num_iter` setting in both Wannier90 setup and convergence validation.

**Architecture:** Thread one positive integer from `main_dft` through the existing Gamma-library setup and run APIs.  It replaces both hard-coded 400 literals without changing the localization algorithm or tolerance.

**Tech Stack:** Fortran 2008, MPI, Wannier90 library, Python route checks, CMake.

---

### Task 1: Add the failing route contract

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the failing test**

Require both production calls to pass `wannier_num_iter`, require the setup
routine to write its integer argument, and require the run routine to pass that
argument to `validate_dg_w90_convergence_log`.  Reject the fixed 400 literals.

**Step 2: Run test to verify it fails**

Run: `python3 tests/dg/check_dg_overlapping_wannier_route.py`

Expected: FAIL because production still contains `num_iter = 400` and validates
against literal 400.

### Task 2: Thread the existing setting through production

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/main_dft.f90`

**Step 1: Implement the minimal API change**

Add a positive `num_iter` input to `setup_dg_w90_gamma_library` and
`run_dg_w90_gamma_library`.  Include it in collective contract validation,
write it to the `.win`, and use it for convergence-log validation.  Import and
pass the existing `wannier_num_iter` from `main_dft`.

**Step 2: Run focused checks**

Run: `python3 tests/dg/check_dg_overlapping_wannier_route.py`

Expected: PASS.

Run: `python3 tests/dg/check_sawf_dmn_format.py --build-dir /Users/otobetoshihito/SALMON-dev/verification/20260818-onepass-production-build`

Expected: PASS.

Run: `git diff --check`

Expected: no output.

### Task 3: Rebuild and prepare Si64 verification

**Files:**
- Modify only the copied verification `inputfile`, not repository fixtures.

**Step 1: Rebuild SALMON serially**

Run: `cmake --build /Users/otobetoshihito/SALMON-dev/verification/20260818-onepass-production-build --target salmon -j1`

Expected: `Built target salmon`.

**Step 2: Prepare a clean run directory**

Copy the known Si64 input and material files, add `wannier_num_iter=1000` to the
existing Wannier namelist, and preserve MPI 8 / OMP 1.

**Step 3: Run and monitor**

Verify that the generated `.win` contains `num_iter = 1000`, Wannier90 reports
convergence before the limit, and SALMON proceeds into post-Wannier processing.
