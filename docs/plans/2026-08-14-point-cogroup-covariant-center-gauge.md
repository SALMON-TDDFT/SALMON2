# Point-Cogroup Covariant Center Gauge Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Rotate each translation-character internal space into a deterministic periodic-position gauge that is transported consistently by point-cogroup and Gamma actions and passes the existing affine center-orbit gate.

**Architecture:** Build and canonicalize three small projected periodic-position matrices in one reference character sector, use the LCFO projected operator only inside position-degenerate blocks, and transport the resulting block gauge through measured point and Gamma intertwiners.  Stream one character orbit at a time and retain only row-owned sector data plus `O(m^2)` dense matrices.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, EigenExa-adjacent row ownership, Python MPI runners and source route contracts, CMake.

---

### Task 1: Add projected periodic-position tuple REDs

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_w90_mpi.py`

1. Import a new `build_dg_sector_periodic_position_tuple` API.
2. Construct a row-owned two-channel sector with known complex periodic phases.
3. Compare all three returned `m x m` matrices with a direct global sum.
4. Rotate the input sector by a complex unitary and require `Z_a -> U^H Z_a U`.
5. Add duplicate ownership, rank-disagreeing metadata, nonfinite, finite-huge, and receipt checks.
6. Print a decomposition-independent tuple fingerprint and compare MPI 1/2/4/8.
7. Run `python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py` and verify RED compilation failure.
8. Commit as `test(dg): cover sector periodic position tuple`.

### Task 2: Implement the distributed periodic-position tuple

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`

1. Add the public API accepting `comm`, row IDs, global row count, row-owned sector, three local unit phases, tolerance, and provenance fingerprint.
2. Establish collective dimension/tolerance/provenance agreement and exactly-once row ownership before indexed work.
3. Bound finite inputs before products, preflight all extents/bytes with checked `int64`, and use collective allocation consensus/cleanup.
4. Compute `Z_a(i,j)=sum_r conjg(C(r,i))*phase_a(r)*C(r,j)` with local BLAS and `MPI_Allreduce` over only `3*m*m` complex values.
5. Measure Gram-I and tuple finite/Hermitian-conjugate contracts, generate a quantized gauge-covariant spectrum/projector fingerprint, and publish persistent/transient workspace receipts.
6. Run W90 MPI 1/2/4/8 and `git diff --check`; commit as `feat(dg): build sector periodic position tuple`.

### Task 3: Add internal canonical-gauge REDs

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`

1. Import `canonicalize_dg_sector_periodic_position_gauge`.
2. Add a two-center tuple whose input columns are mixed by a known complex unitary; require direct equality of canonical frames after input rotation.
3. Add a position-degenerate block resolved by a nondegenerate LCFO operator.
4. Add an exact multiplet with both discriminators degenerate; require identical final projector/fingerprint under arbitrary internal rotation rather than arbitrary column equality.
5. Add inconsistent operator covariance, split-cluster, rank disagreement, and finite-huge REDs.
6. Run W90 MPI fixture and verify RED failure before implementation.

### Task 4: Implement deterministic internal canonicalization

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Form a bounded Hermitian discriminator from fixed irrational coefficients multiplying Hermitian real/imaginary parts of the three position matrices.
2. Diagonalize with checked `ZHEEV`, collectively agree every adjacent spectral cluster boundary, and measure eigensystem residuals.
3. Within each degenerate block, project/diagonalize the LCFO operator and collectively agree its clusters.
4. Phase-fix nondegenerate columns using the existing smallest-global-row maximum-amplitude convention.
5. Preserve unresolved exact multiplets as common projectors and produce deterministic block metadata; do not split them with LAPACK column order.
6. Apply the `m x m` unitary to row-owned sector rows with BLAS, validate Gram-I and tuple covariance, and bind all inputs/blocks to the output fingerprint.
7. Run W90 MPI 1/2/4/8; commit as `feat(dg): canonicalize periodic position gauge`.

### Task 5: Add point-character transport and Gamma REDs

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Build a Z4 fixture where a point representative maps a character to its non-self-inverse partner with a known internal action.
2. Require the transported position tuple to satisfy integer-rotation covariance with the correct left/right convention.
3. Build a Z2 x Z2 fixture with self-conjugate sectors and require the existing Gamma fixed-frame residual.
4. Corrupt the character permutation, point action orientation, Gamma pairing, and exact-multiplet metadata independently and require collective rejection.
5. Require inverse-transformed values to be Gamma-real and the center-gauge diagnostic leakage to vanish.
6. Run both focused MPI runners and verify RED failures.

### Task 6: Implement one-orbit covariant gauge transport

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

1. Add a compact point-character orbit descriptor derived from canonical character indices and point-cogroup action.
2. Canonicalize one reference sector, then transport its internal block gauge with the measured point intertwiner to every character in that orbit.
3. Process every conjugate pair exactly once; generate non-self-conjugate partners through Gamma and use the fixed-frame branch for self-conjugate sectors.
4. Validate periodic-position covariance on point generators, Gamma pairing, orbit completeness, and duplicate/missing scheduling before returning any sector.
5. Keep only one aligned sector/conjugate pair and `O(m^2)` matrices live; publish checked collective peak receipts and aggregate provenance.
6. Run W90 and construction MPI 1/2/4/8; commit as `feat(dg): transport covariant point center gauge`.

### Task 7: Connect production with route RED first

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `src/gs/main_dft.f90`

1. Add a route RED requiring position tuple construction and point-character transport after reference-sector anchoring and before inverse accumulation.
2. Require one-character-orbit scheduling, Gamma pairing before accumulator insertion, and no all-sector/full-rank replicated arrays.
3. Run the route checker and observe failure.
4. Connect the production loop using existing local periodic phases, LCFO operator, character catalog, point representatives/cocycle, and Gamma receipts.
5. Aggregate pre/post spread, rotation norm, position covariance, Gamma, center leakage, workspace, and provenance receipts across all orbits.
6. Keep the existing post-inverse factored cocycle and center gates authoritative.
7. Run route checker GREEN and commit as `feat(dg): apply covariant point center gauge in production`.

### Task 8: Focused verification and review

**Files:**
- Verify only.

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_obsolete_dg_routes_removed.py
cmake --build /Users/otobetoshihito/SALMON-dev/verification/20260814-core-ownership-fix-build -j 1
git diff --check
```

Expected: all PASS.  Review collective branch agreement, default/MPI count bounds, allocation cleanup, action direction, Gamma scheduling, exact multiplets, and one-sector memory before a production run.

### Task 9: Si64 validation

1. Create a fresh verification directory and record commit/binary hash.
2. Run MPI 8 with `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`.
3. Monitor CPU, per-rank RSS, and log progress without changing rank count.
4. Require point-center monomial/leakage diagnostics, affine center closure, factored cocycle proof, Gamma reality, and spread receipts to pass.
5. If the gate still fails, stop at the diagnostic evidence; do not relax tolerance or stack another repair without revisiting the architecture.
