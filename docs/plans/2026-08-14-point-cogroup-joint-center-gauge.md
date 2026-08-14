# Point-Cogroup Joint Center Gauge Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the noncovariant fixed-direction periodic-position gauge with a point-cogroup-invariant joint center gauge that passes the unchanged affine center-orbit gate.

**Architecture:** Jointly diagonalize the six Hermitian real/imaginary components of the weighted periodic-position tuple using bounded small-matrix Jacobi sweeps. Resolve repeated-center blocks with the LCFO operator, preserve exact multiplets, and transport the resulting block gauge through existing point-character and Gamma intertwiners one orbit at a time.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, checked `int64` receipts, Python route contracts, CMake.

---

### Task 1: Add joint-center REDs

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_w90_mpi.py`

1. Import a new `jointly_canonicalize_dg_sector_periodic_position_gauge` API.
2. Construct two distinct periodic centers, rotate their columns by a complex unitary, and require the same sorted centers and direct canonical frame.
3. Rotate the Cartesian tuple with a nontrivial integer point rotation and require the same center set/projector fingerprint.
4. Add a coincident-center block resolved by LCFO and an exact LCFO-degenerate multiplet checked only by projector.
5. Add nonconvergence, rank-disagreement, duplicate ownership, nonfinite, finite-huge, and unsafe-workspace REDs.
6. Print and compare the joint-center fingerprint on MPI 1/2/4/8.
7. Run `python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py`; expect missing-symbol compilation failure.
8. Commit `test(dg): cover joint periodic center gauge`.

### Task 2: Implement bounded joint diagonalization

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Add the public API with row ownership, weighted tuple, LCFO operator, tolerance, and provenance inputs.
2. Establish collective metadata/payload agreement and preallocation checked receipts.
3. Form six Hermitian matrices and implement exact complex 2x2 Jacobi updates with no division inside the row/column loops.
4. Bound sweeps, collectively agree convergence, and publish initial/final objective and maximum update.
5. Sort columns by the three periodic center phases; cluster adjacent periodic distances.
6. Diagonalize LCFO inside repeated-center blocks; preserve exact multiplet projectors.
7. Phase-fix nondegenerate columns, validate Gram and tuple covariance, and stream an `O(N*m)` fingerprint.
8. Run W90 MPI 1/2/4/8 and `git diff --check`; commit `feat(dg): jointly canonicalize periodic centers`.

### Task 3: Add point-character transport REDs

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Add a directed Z4 point action mapping a character to its non-self-inverse partner.
2. Require transported tuples and sorted centers to obey the supplied integer Cartesian rotation.
3. Add Z2xZ2 self-conjugate Gamma-real coverage and a non-self-conjugate Gamma pair.
4. Corrupt the action direction, character permutation, exact-multiplet metadata, and Gamma pairing independently.
5. Require zero center-block leakage after inverse reconstruction.
6. Run both MPI runners and observe RED failures; commit `test(dg): cover point covariant center transport`.

### Task 4: Implement one-orbit center transport

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

1. Build the compact point-character orbit from canonical catalog metadata.
2. Jointly localize one reference sector and transport its block gauge with measured point intertwiners.
3. Schedule every character/conjugate exactly once; use Gamma generation/fixed-real branches.
4. Validate position covariance on point generators, orbit completeness, duplicate/missing scheduling, and final Gamma residual.
5. Return one row-owned sector/pair at a time with checked `O(N*m+m^2)` peak receipts.
6. Run W90 and construction MPI 1/2/4/8; commit `feat(dg): transport joint center gauge`.

### Task 5: Replace the production fixed-direction gauge

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `src/gs/main_dft.f90`

1. Add a route RED forbidding `canonicalize_dg_sector_periodic_position_gauge` in production and requiring the joint/transport APIs before inverse accumulation.
2. Run the route checker and observe failure.
3. Remove the fixed-direction production call while retaining the weighted tuple builder.
4. Connect joint reference localization and one-orbit point/Gamma transport.
5. Aggregate joint objective, covariance, Gamma, workspace, and provenance receipts.
6. Keep factored cocycle and affine center gates authoritative.
7. Run route checker and focused MPI tests; commit `feat(dg): apply joint point center gauge in production`.

### Task 6: Verify and run Si64

**Files:**
- Verify only.

1. Run W90 MPI 1/2/4/8.
2. Run construction MPI 1/2/4/8.
3. Run route and obsolete-route contracts.
4. Build Release with `-j 1` and run `git diff --check`.
5. Create a fresh Si64 directory recording commit and binary hashes.
6. Run MPI 8 with `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and unchanged tolerance.
7. Require the joint canonicalization to complete in bounded time, factored closure to pass, and operation-2 center residual/leakage to fall below tolerance.
8. If it fails, stop at the measured joint objective/covariance evidence; do not relax tolerance.
