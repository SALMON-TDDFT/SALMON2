# LCFO Boundary-Calibrated Occupied Symmetry Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Retain the coherent rank-128 LCFO occupied space while distinguishing fragment-boundary stitching error from genuine interior symmetry breaking.

**Architecture:** Measure the full-group rank-128 projection residual on unique cores, split it into stencil-width boundary and interior components, and calibrate the boundary allowance from LCFO face value/gradient mismatch.  If the physical gates pass, polar-correct and group-synchronize the measured representation inside rank 128; never promote stitching residuals to occupied states.

**Tech Stack:** Fortran 2008, MPI, LAPACK, EigenExa, Python source/evidence contracts.

---

### Task 1: Classify LCFO covariance residuals

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Add a RED MPI case with the same perturbation placed first on a boundary point and then on an interior point.
2. Verify RED on 1/2/4/8 ranks.
3. Add a distributed rank-fixed projection-residual routine returning total, boundary, and interior norms per operation.
4. Build the boundary mask from unique-core indices and the derivative-stencil radius; do not remove boundary points from quadrature.
5. Run the fixture on 1/2/4/8 ranks and `git diff --check`.

### Task 2: Calibrate and gate the residual

**Files:**
- Modify: `src/gs/dc/lcfo.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Add RED cases requiring the allowance to scale with measured face value/gradient mismatch and rejecting an equal interior residual.
2. Return LCFO face value and gradient mismatch from the in-memory occupied contribution path.
3. Propagate the measured mismatch through grid spacing and stencil radius to a dimensionless boundary leakage allowance.
4. Gate strict interior leakage, calibrated boundary leakage, density covariance, and inversion-odd density separately.
5. Report every raw diagnostic and run focused construction/route tests.

### Task 3: Correct within rank 128

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Add RED cases for nonunitary measured operations, product-table drift, projector change, and a cut symmetry block.
2. Metric-orthonormalize the LCFO occupied rows without changing their projector.
3. Polar-project the measured rank-128 operations and synchronize them to the full product table.
4. Gate correction norm, projector change, density change, and occupied energy change against the measured boundary baseline.
5. Add optional projection seeds only after the rank-128 gate and select complete blocks at rank 384.
6. Run construction, symmetry, metric, projection, and route fixtures on 1/2/4/8 ranks.

### Task 4: Genuine evidence and reviews

**Files:**
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`

1. Run strict undisplaced ideal-Si64 DC-SCF, LCFO+EigenExa, and report the calibrated diagnostics.
2. Complete specification review and code-quality review; resolve every Critical/Important finding.
3. Run rank-384 Wannier localization and publish V3 only after all gates pass.
4. Run field-off, linear-response, and long-pulse Exp RT; generate polarization-derived semi-log HHG with even peak/dip/slope classifications.
5. Perform a clean-first committed-HEAD EigenExa build, final focused verification, and repeat both reviews.
6. Commit evidence and push the identical head to origin and upstream.

