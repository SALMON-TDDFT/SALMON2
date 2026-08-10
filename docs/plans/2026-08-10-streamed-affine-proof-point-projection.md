# Streamed Affine Proof and Point-Group Projection Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove the replicated full-affine representation from the overlapping-Wannier production route while preserving a streamed proof of every affine operation and dense projection under only the exact full-system fixed-center point subgroup.

**Architecture:** Treat the complete affine space group as a streamed validation layer and the maximal common-fixed-center point subgroup as a separate bounded dense operator-projection layer.  Orthonormalize the already affine-closed LCFO seed space without regenerating every orbit, run and canonicalize Wannier90, validate all affine operations through row-owned streaming, and publish both proof layers in V3.

**Tech Stack:** Fortran 2008, MPI, EigenExa, ScaLAPACK, spglib, Wannier90 3.1 library mode, Python source contracts, CMake clean overlays.

---

## Mandatory execution discipline

For every task:

1. use `@test-driven-development` and capture a genuine RED before production edits;
2. run focused tests on MPI 1/2/4/8 where applicable;
3. perform specification and code-quality reviews;
4. resolve all Critical and Important findings and repeat affected tests;
5. build a clean `git archive HEAD` overlay with only the task diff and committed parent prerequisites;
6. keep normal DC LCFO+EigenExa and generalized-eigenvalue Exp-only V3 RT unchanged; and
7. do not push until the complete parent plan has passed genuine Si64 GS/LR/HHG acceptance.

### Task AP1: Orthonormalize an already affine-closed distributed seed space

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the failing tests**

Add a fixture with a seed space whose measured full-group residual is below tolerance.  Require a new distributed orthonormalization routine to return exactly the seed rank, preserve its projector, and use no symmetry-image orbit expansion.  Add malformed metric, rank-inconsistent seed count, and non-finite input rejection.  Extend the route contract to forbid both `Nsym**2*Ngrid` action revalidation and a production call that regenerates the full orbit after the affine residual has passed.

**Step 2: Run RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because production still calls `build_dg_distributed_symmetry_closed_basis` with all affine operations and the new orthonormalization API is absent.

**Step 3: Implement the minimal distributed path**

Use the existing direct cyclic EigenExa metric machinery to compute the seed metric and inverse square root without a replicated dense input.  Apply the inverse square root in orbital tiles to produce an orthonormal row-distributed/core-distributed basis.  Accept this shortcut only after the previously computed maximum full-affine LCFO residual is within `dg_ow_symmetry_tolerance`; otherwise reject rather than falling back to orbit expansion.

Delete the redundant pairwise map-composition/grid proof from the generic builder.  Retain collective permutation validation and product-table range validation for its fixtures and remaining nonproduction callers.

**Step 4: Focused verification**

Run construction and metric fixtures on MPI 1/2/4/8 and the route checker.  Compare projectors with the existing small dense reference and verify bounded workspace receipts.

**Step 5: Reviews and commit**

Review preservation of occupied and complete-s+p ranks, Gamma-real handling, rank-count determinism, and allocation overflow.  Resolve all Critical/Important findings, perform the clean overlay build, then commit only AP1 files.

### Task AP2: Stream the full affine proof and materialize only the fixed-center point subgroup

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_symmetry_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the failing tests**

For a small affine group with nontrivial translations, compare streamed identity/unitarity/closure/fingerprint receipts with a dense reference.  Separately construct the maximal operations sharing a supplied full-system center, require a closed point subgroup and exact inversion, and compare scalar/vector operator projection with the dense reference.  Add a source RED forbidding `Nsym*Nwann**2` allocation and a `global_retained_representation(:,:,Nsym)` declaration in production.

**Step 2: Run RED**

Run construction and symmetry fixtures on MPI 1/2/4/8 plus the route checker.  Expected: FAIL because `main_dft` still materializes `global_symmetry_overlap`, `global_candidate_raw`, and `global_retained_representation` for all affine operations.

**Step 3: Implement streamed proof**

Retain the full affine point maps, multiplication/factorization metadata, and row-owned post-MLWF validator.  Stream one operation and orbital row tile at a time, update identity/unitarity/closure maxima and a process-count-independent fingerprint, then release the tile.  Validate closure through the factored translation/point-cogroup generator relations rather than all `Nsym**2` products over the grid.

**Step 4: Implement fixed-center subgroup projection**

Derive the center from the full atomic catalog, not fragments.  Select every exact operation fixing that center modulo a lattice vector, deduplicate rotations, build its closed product table, and materialize at most 48 dense representation matrices.  Use this representation for overlap/Hamiltonian/position/velocity projection; shift position by the same center.  Keep translation and cocycle receipts in the streamed proof rather than the dense tensor.

**Step 5: Focused verification**

Run symmetry, construction, operator, solver, SCF, and coefficient-RT fixtures on MPI 1/2/4/8.  Require dense-reference equivalence for small groups, exact inversion retention, process-count-independent fingerprints, and no forbidden allocation.

**Step 6: Reviews and commit**

Review group composition order, affine pullback convention, off-fragment centers, Gamma phase, position-origin handling, workspace accounting, and error paths.  Resolve all Critical/Important findings, perform a clean full-feature overlay build, and commit only AP2 files.

### Task AP3: Bind both proof layers to V3 and rerun genuine ideal Si64

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`

**Step 1: Write publication RED tests**

Require V3 to contain a full-affine streamed fingerprint, affine identity/unitarity/closure receipts, translation-subgroup order, point-cogroup order, fixed-center subgroup fingerprint/order/center, inversion receipt, and separate workspace peaks.  Reject zero, nonfinite, rank-inconsistent, over-tolerance, over-limit, or missing-inversion evidence for the centrosymmetric fixture.

**Step 2: Run RED**

Run the checkpoint fixture on MPI 1/2/4/8 and the route checker.  Expected: FAIL on absent two-layer provenance.

**Step 3: Implement V3 serialization and gates**

Extend manifest write/read, expected size, digest, broadcast, replicated-payload consistency, rejection code, and production population.  Keep the existing MLWF backend/version, M/A and transform fingerprints, spreads, coordinator bytes, canonical flag, and generalized-eigenvalue Exp reader boundary.

**Step 4: Clean overlay and genuine Si64 GS**

Build committed HEAD plus reviewed AP1--AP3/W4 diffs with MPI, ScaLAPACK, EigenExa, spglib, and Wannier90 enabled.  Build EigenExa with `-j1`, then SALMON `--clean-first -j4`.  Run ideal undisplaced Si64 on 8 MPI ranks and 1 OpenMP thread.  Require 384 retained/128 occupied states, full affine proof, fixed-center inversion, bounded memory, canonical MLWF provenance, reconstructed-GS receipts, and accepted V3.

**Step 5: LR and HHG acceptance**

Run field-off, impulse x/y/z and half-amplitude, weak laser, 10-cycle strong laser x/y/z, and half-dt convergence cases using only generalized-eigenvalue Exp coefficient RT.  Compute LR/HHG spectra from polarization; use current only for `dP/dt`.  Produce the absolute-path semilog figure and classify H2/H4 as peak/dip/slope, requiring even-harmonic suppression relative to H3 for ideal Si64.

**Step 6: Final reviews and commit**

Run every retained focused fixture on MPI 1/2/4/8, route/removal contracts, genuine Si64 checkers, `git diff --check`, and a final committed-HEAD clean overlay.  Resolve all Critical/Important findings before committing AP3/W4.  Continue with parent Tasks 5--6; push to both `origin` and `upstream` only after the entire plan passes.
