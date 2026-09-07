# Wannier Center Action Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Distinguish a genuine Wannier90 symmetry failure from an affine-action convention or origin mismatch when the Si64 Wannier-center orbit gate fails.

**Architecture:** Keep the existing forward-action acceptance gate unchanged. On failure only, compute diagnostic nearest-neighbour residuals for the inverse affine action and for an origin-offset fit, then include those values and the failing source/mapped coordinates in the message. No production state is changed.

**Tech Stack:** Fortran 2008, existing DG construction MPI fixture, Python MPI runner.

---

### Task 1: Lock the diagnostic contract with a RED

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Extend the existing broken-center-orbit assertion to require forward and inverse residual labels and source/mapped coordinates.
2. Run `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py` and confirm failure because the labels are absent.

### Task 2: Add failure-only action diagnostics

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

1. On the first failed matching operation, compute the inverse affine image using `R^{-1}(x-tau)` for unimodular integer rotations.
2. Compute its nearest periodic target residual and report both residuals plus source and forward-mapped coordinates.
3. Keep the forward matching result as the sole accept/reject criterion.
4. Run the focused construction MPI test on 1/2/4/8 ranks.

### Task 3: Reproduce Si64 and classify the failure

**Files:**
- No source changes.

1. Rebuild the existing lightweight executable.
2. Rerun the Si64 case to the immediate post-Wannier90 center gate.
3. Compare forward and inverse residuals. If neither is small, inspect a common origin offset before changing symmetry logic.
4. Run `git diff --check` and the route checker.
