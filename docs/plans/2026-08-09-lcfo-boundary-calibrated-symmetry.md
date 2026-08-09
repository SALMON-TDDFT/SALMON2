# LCFO Boundary-Calibrated Occupied Symmetry Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Generate all 384 Wannier states from one coherent LCFO/EigenExa coefficient space while distinguishing fragment-boundary stitching error from genuine interior symmetry breaking.

**Architecture:** Ask LCFO/EigenExa for 384 coefficient columns, evaluate the corresponding fragment-basis contributions only on core plus buffer, and sum them onto unique cores by physical grid ID.  Measure full-group residuals inside this fixed rank-384 space, calibrate boundary error independently, and symmetry-correct only within rank 384.  The first 128 columns define the occupied projector; fragment-local complement states are forbidden.

**Tech Stack:** Fortran 2008, MPI, LAPACK, EigenExa, Python source/evidence contracts.

---

### Task 1: Classify LCFO covariance residuals

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Add RED source and numerical cases requiring LCFO to return more coefficient columns than `dc%nstate_tot`, up to the requested rank 384.
2. Verify RED before changing LCFO allocation or EigenExa extraction bounds.
3. Generalize the in-memory LCFO contribution output to the requested retained count while leaving normal DC LCFO output unchanged.
4. Remove the fragment-local occupied/complement construction from the OW production route; use the first 384 LCFO rows and the first 128 as occupied.
5. Add a RED MPI case with the same perturbation placed first on a boundary point and then on an interior point.
6. Add a distributed rank-fixed projection-residual routine and stencil-width boundary mask.
7. Run the fixture on 1/2/4/8 ranks, the route contract, full build, and `git diff --check`.

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

### Task 3: Correct within rank 384

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Add RED cases for nonunitary measured operations, product-table drift, projector change, and a cut symmetry block.
2. Metric-orthonormalize the LCFO rank-384 rows while verifying inclusion of the first 128 occupied rows.
3. Polar-project the measured rank-384 operations and synchronize them to the full product table.
4. Gate correction norm, projector change, density change, and occupied energy change against the measured boundary baseline.
5. Use atomic projection overlap only as a symmetry-averaged localizer inside rank 384; never add projector or fragment-local rows.
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
