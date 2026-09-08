# Wannier90 Symmetry-Projected Line Search Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make Wannier90's line-search derivative describe the symmetry-projected direction that is actually applied.

**Architecture:** Extend the existing bundled-Wannier90 source patch at the projection boundary.  Preserve the upstream optimizer and add only post-projection derivative/norm recomputation plus a safe CG reset when projection makes the direction non-descending.

**Tech Stack:** CMake source patching, Fortran 2008, MPI reductions, Python route tests, bundled Wannier90 3.1.

---

### Task 1: Lock the source-order contract

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Test: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Add assertions requiring projected slope recomputation after `sitesym_symmetrize_gradient` and before the fixed-step/line-search branch.
2. Require a projected-direction diagnostic and CG-reset guard.
3. Run the route test and verify it fails because the patch does not yet contain those operations.

### Task 2: Implement the minimal post-projection correction

**Files:**
- Modify: `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`
- Test: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Extend the exact source replacement after the symmetry projection.
2. Recompute `doda0 = -Re<gradient,direction>/(4*wbtot)` using owned slices and an MPI all-reduction.
3. Compute/log the projected direction norm.
4. If the projected direction is uphill, reset CG, project the steepest-descent direction, and recompute once.
5. Run the route test and verify it passes.

### Task 3: Build and inspect generated Wannier90

**Files:**
- Inspect generated: `wannier90/src/wannier90-project/src/wannierise.F90`

1. Reconfigure/rebuild the existing production build.
2. Verify the generated source order and build success.
3. Run focused SALMON/Wannier90 route tests.

### Task 4: Short same-checkpoint A/B validation

**Files:**
- Create: a new directory under `/Users/otobetoshihito/SALMON-dev/verification/`

1. Copy the established Si64 checkpoint/replay inputs without changing the physical system.
2. Enable verbose line-search diagnostics and run a bounded diagnostic segment, not a production convergence run.
3. Compare projected slope, accepted spread trajectory, gradient diagnostics, runtime, and memory with the previous trajectory.
4. Continue to a long run only if the diagnostic confirms direction/derivative consistency and improves stability.

