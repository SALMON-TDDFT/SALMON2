# Wannier90 Local Gradient Synchronization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Ensure Wannier90's symmetry-projected search direction is the direction actually used by its distributed line search.

**Architecture:** Extend the existing external-project patch so the distributed `cdq_loc` direction is gathered into `cdq`, projected, and copied back to each rank's local slice. Protect the data-flow ordering with a source-level regression, then verify the rebuilt library and Si64 runtime gate.

**Tech Stack:** CMake script-mode patching, Fortran 2008, MPI, Wannier90, Python regression checks.

---

### Task 1: Add a failing patch data-flow regression

**Files:**
- Modify: `tests/dg/check_sawf_dmn_format.py`
- Test: `tests/dg/check_sawf_dmn_format.py`

**Step 1:** Add assertions that the patch contains `comms_gatherv`,
`comms_bcast`, `sitesym_symmetrize_gradient(2, cdq)`, and the owned-slice copy
to `cdq_loc`, in that order.

**Step 2:** Run `python3 tests/dg/check_sawf_dmn_format.py` and verify it fails
because the synchronization patch is absent.

### Task 2: Patch the Wannier90 optimizer

**Files:**
- Modify: `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`

**Step 1:** Extend the CMake patch to replace the lone gradient projection call
with gather, broadcast, projection, and local-slice restoration.

**Step 2:** Run `python3 tests/dg/check_sawf_dmn_format.py` and verify it passes.

**Step 3:** Apply the patch to a clean Wannier90 source copy and inspect the
resulting `wannierise.F90` data flow.

### Task 3: Rebuild and run focused verification

**Files:**
- No new source files.

**Step 1:** Rebuild the external Wannier90 library and SALMON overlay.

**Step 2:** Run W90 MPI tests on 1/2/4/8 ranks.

**Step 3:** Run construction MPI tests on 1/2/4/8 ranks and the route checker.

**Step 4:** Run `git diff --check`.

### Task 4: Validate Si64 behavior

**Files:**
- Runtime output only under the external verification directory.

**Step 1:** Rerun Si64 with the same MPI rank count and input.

**Step 2:** Confirm Wannier90 converges and the post-Wannier generator
covariance gate passes.  If it fails, retain the checkpoint and report the new
per-generator defects before making another change.

