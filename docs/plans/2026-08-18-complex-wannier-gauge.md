# Complex Wannier Gauge Acceptance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Accept physically equivalent complex unitary Wannier90 gauges.

**Architecture:** Remove element-wise-real gates from the immediate result validator and production applicator. Canonicalize only arbitrary output-column phases through a deterministic spatial pivot, applying the same phase to the transform, values, and gradients. Preserve all gauge-invariant numerical validation and downstream symmetry checks.

**Tech Stack:** Fortran 2008, MPI fixture, Python route checker.

---

### Task 1: Add a complex-unitary RED

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Validate `diag(1,i)` as a finite unitary transformation and require success.
2. Require the obsolete Gamma-real error path to be absent.
3. Run the focused route and MPI tests and observe failure under the old gate.

### Task 2: Remove only the gauge-dependent rejection

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Delete the imaginary-part defect and its scale calculation.
2. Keep finite, spread, unitarity, and spread-improvement gates unchanged.
3. Run the route checker and MPI fixture on 1/2/4/8 ranks.

### Task 3: Verify production integration

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Add a production-applicator fixture with a genuinely complex unitary transform and nonzero complex gradients.
2. Verify that deterministic pivot phases are applied identically to transform columns, values, and gradients.
3. Remove the applicator's imaginary-transform and real-pivot/sign-only assumptions; use `conjg(pivot)/abs(pivot)` instead.
4. Run the focused MPI test on 1/2/4/8 ranks.

### Task 4: Verify and review

1. Run the SAWF/DMN format checker.
2. Incrementally rebuild SALMON.
3. Run `git diff --check` and commit only the scoped files.
4. Request independent code review of the full production path.
5. Re-run the Si64 case only after Critical/Important review findings are resolved.
