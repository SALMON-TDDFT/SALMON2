# Exact Buffered-Fragment Crystallographic Symmetry Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Detect and enforce only the exact instantaneous crystallographic subgroup of each buffered fragment so arbitrary nonmagnetic materials retain correct tensor covariance without erasing zero-point or thermal symmetry breaking.

**Architecture:** Reuse spglib as the crystallographic operation oracle, normalize its operations into one audited catalog, and filter that catalog independently for each core-plus-buffer fragment. Build metric-unitary retained-basis representations and apply local Reynolds projection to scalar and polar-vector matrices before publishing the unchanged V3 checkpoint; C1 is a valid no-projection result.

**Tech Stack:** Fortran 2008, C spglib wrapper/API, MPI, ScaLAPACK, EigenExa, Python 3 contract and production runners, CMake/CTest.

---

### Task 1: Freeze the failed physical matrix and exact-symmetry policy

**Files:**
- Add: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`
- Add: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Add: `tests/dg/data/si64_overlapping_wannier_rt/input_gs.in`
- Add: `tests/dg/data/si64_overlapping_wannier_rt/atom.dat`
- Add: `tests/dg/data/si64_overlapping_wannier_rt/input_fieldoff.in`
- Add: `tests/dg/data/si64_overlapping_wannier_rt/input_impulse.in`
- Add: `tests/dg/data/si64_overlapping_wannier_rt/input_laser_weak.in`
- Add: `tests/dg/data/si64_overlapping_wannier_rt/input_laser_hhg.in`
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`

**Step 1: Preserve RED evidence**

Record the genuine Si64 checkpoint provenance and the completed eleven-case
matrix.  Record that field-off drift, impulse linearity, `J` versus `dP/dt`,
and time-step convergence pass while cubic-axis and odd/even HHG gates fail.

**Step 2: Run the focused checker and confirm RED**

Run:

```bash
python3 tests/dg/check_si64_overlapping_wannier_response_hhg.py \
  /tmp/si64-response-hhg-parent.3bSpnu/matrix
```

Expected: FAIL at cubic-axis equivalence or odd-harmonic dominance.  Do not
weaken the physical tolerances.

**Step 3: Verify fixtures and source contracts**

Run:

```bash
python3 -m py_compile tests/dg/check_si64_overlapping_wannier_response_hhg.py \
  tests/dg/run_si64_overlapping_wannier_response_hhg.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_obsolete_dg_routes_removed.py
git diff --check
```

Expected: PASS except the intentional production physics RED from Step 2.

**Step 4: Review and commit**

Perform specification and code-quality reviews.  Resolve every Critical and
Important finding, rerun affected checks, then commit:

```bash
git add tests/dg/check_si64_overlapping_wannier_response_hhg.py \
  tests/dg/run_si64_overlapping_wannier_response_hhg.py \
  tests/dg/data/si64_overlapping_wannier_rt \
  docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md
git commit -m "test(dg): expose Si64 response symmetry failure"
```

### Task 2: Survey the crystallographic operation catalog

**Files:**
- Modify: `src/gs/dc/lcfo_wannier_spglib.c`
- Modify: `src/gs/dc/lcfo_wannier_sawf.f90`
- Modify: `tests/dg/fixtures/sawf_spglib/driver.F90`
- Modify: `tests/dg/fixtures/sawf_spglib/mock_spglib.c`
- Add: `tests/dg/check_dg_crystallographic_point_groups.py`
- Add: `tests/dg/fixtures/crystallographic_point_groups.json`

**Step 1: Write catalog RED tests**

Add fixtures for all 32 crystallographic point groups with integer fractional
rotations, determinants, operation orders, and multiplication tables.  Include
representative fractional translations from nonsymmorphic space groups.  Test
that the spglib wrapper returns the dataset identity, point-group symbol, and
normalized operations without assuming orthogonal Cartesian axes.

**Step 2: Run RED**

Run:

```bash
python3 tests/dg/check_dg_crystallographic_point_groups.py
python3 tests/dg/check_sawf_spglib_backend.py
```

Expected: FAIL because the catalog metadata and complete survey API do not yet
exist.

**Step 3: Implement the minimal catalog API**

Extend the existing wrapper rather than adding a second symmetry backend.
Normalize rotations and translations through `lcfo_wannier_sawf.f90`; expose
stable point-group metadata and deterministic ordering.  Validate identity,
duplicates, inverse, determinant, and closure.

**Step 4: Focused verification**

Rerun both commands from Step 2 and `git diff --check`.  Expected: PASS.

**Step 5: Reviews and commit**

Complete specification and code-quality reviews, resolve all Critical and
Important findings, rerun affected checks, then commit:

```bash
git add src/gs/dc/lcfo_wannier_spglib.c src/gs/dc/lcfo_wannier_sawf.f90 \
  tests/dg/fixtures/sawf_spglib tests/dg/check_dg_crystallographic_point_groups.py \
  tests/dg/fixtures/crystallographic_point_groups.json
git commit -m "feat(dg): survey crystallographic point groups"
```

### Task 3: Select the exact instantaneous subgroup per buffered fragment

**Files:**
- Modify: `src/gs/dc/lcfo_wannier_sawf_band.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Add: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Add: `tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90`
- Add: `tests/dg/run_dg_overlapping_wannier_fragment_symmetry_mpi.py`

**Step 1: Write subgroup RED tests**

Cover exact triclinic, monoclinic, orthorhombic, tetragonal, trigonal,
hexagonal, and cubic fragments; multiple species; improper rotations;
periodic images; fragment boundary incompatibility; and different groups in
neighboring fragments.  A roundoff-scale coordinate perturbation must retain
the operation, while a fixed physical displacement must remove it and may
produce C1.  Verify deterministic largest-closed-subgroup selection.

**Step 2: Run RED on 1, 2, 4, and 8 ranks**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_fragment_symmetry_mpi.py
```

Expected: FAIL because the per-fragment exact subgroup API is absent.

**Step 3: Implement exact filtering**

Filter operations using atom/species, core-plus-buffer, grid, and center maps.
Keep atom, boundary, grid, and center residuals separate.  Derive scale-aware
numerical tolerances from machine precision and coordinate conditioning; never
inflate them to recover parent symmetry.  Select a deterministic closed
subgroup containing identity and accept C1 without warning.

**Step 4: Focused verification**

Run the MPI runner, existing SAWF fragment-map tests, route contract, and
`git diff --check`.  Expected: PASS on 1/2/4/8 ranks.

**Step 5: Reviews and commit**

Complete both reviews, resolve all Critical and Important findings, rerun
affected checks, then commit:

```bash
git add src/gs/dc/lcfo_wannier_sawf_band.f90 src/gs/dc/lcfo_flux.f90 \
  src/gs/dc/dg_overlapping_wannier_types.f90 \
  src/gs/dc/dg_overlapping_wannier_symmetry.f90 src/gs/dc/CMakeLists.txt \
  tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_fragment_symmetry_mpi.py
git commit -m "feat(dg): detect exact buffered fragment symmetry"
```

### Task 4: Build and validate retained-basis group representations

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_metric.f90`
- Add: `tests/dg/test_dg_overlapping_wannier_point_group_mpi.f90`
- Add: `tests/dg/run_dg_overlapping_wannier_point_group_mpi.py`

**Step 1: Write representation RED tests**

Use nonorthogonal retained bases and groups containing rotations, mirrors, and
inversion.  Check metric unitarity, multiplication-table closure, inverse
consistency, polar correction, deterministic rank independence, and C1.

**Step 2: Run RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_point_group_mpi.py
```

Expected: FAIL because retained-basis point-group representations are absent.

**Step 3: Implement the minimal representation builder**

Assemble buffered real-space overlaps, apply polar correction, and reject
representations exceeding pre- or post-correction residual gates.  Do not
materialize full-system wavefunctions or introduce a global point group.

**Step 4: Focused verification and reviews**

Run the new runner plus existing construction, metric, and row-storage runners
on 1/2/4/8 ranks.  Perform both reviews, resolve all Critical and Important
findings, rerun affected checks, and commit:

```bash
git add src/gs/dc/dg_overlapping_wannier_symmetry.f90 \
  src/gs/dc/dg_overlapping_wannier_construction.f90 \
  src/gs/dc/dg_overlapping_wannier_metric.f90 \
  tests/dg/test_dg_overlapping_wannier_point_group_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_point_group_mpi.py
git commit -m "feat(dg): represent fragment point groups"
```

### Task 5: Project scalar and polar-vector operators covariantly

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_operators.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_observables.f90`
- Modify: `src/gs/main_dft.f90`
- Add: `tests/dg/test_dg_overlapping_wannier_symmetry_projection_mpi.f90`
- Add: `tests/dg/run_dg_overlapping_wannier_symmetry_projection_mpi.py`

**Step 1: Write covariance RED tests**

Construct noisy scalar `S/H` and polar-vector `X/V` matrices for representative
proper and improper point groups.  Require Reynolds projection to restore
scalar invariance and vector covariance, including inversion signs.  Require
C1 to be byte-preserving and large pre-projection defects to fail.

**Step 2: Run RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_symmetry_projection_mpi.py
```

Expected: FAIL because only Hermiticity gates exist.

**Step 3: Implement local Reynolds projection**

Project each fragment before distributed global assembly, using its own exact
group and retained-basis representation.  Validate pre/post residuals and use
the polar-vector Cartesian convention for `X/V` under improper operations.

**Step 4: Focused verification and reviews**

Run the new runner and all operator, observable, physical-matrix, construction,
and route tests on 1/2/4/8 ranks.  Perform both reviews, resolve all Critical
and Important findings, rerun affected checks, and commit:

```bash
git add src/gs/dc/dg_overlapping_wannier_symmetry.f90 \
  src/gs/dc/dg_overlapping_wannier_operators.f90 \
  src/gs/dc/dg_overlapping_wannier_observables.f90 src/gs/main_dft.f90 \
  tests/dg/test_dg_overlapping_wannier_symmetry_projection_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_symmetry_projection_mpi.py
git commit -m "feat(dg): project fragment-covariant operators"
```

### Task 6: Bind exact symmetry evidence to the V3 checkpoint

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Modify: `src/rt/dg/rt_dg_overlapping_wannier.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`
- Modify: `tests/dg/run_rt_dg_overlapping_wannier_mpi.py`

**Step 1: Write checkpoint RED tests**

Require V3 magic to remain unchanged while the manifest/fingerprint binds
per-fragment operations, multiplication tables, tolerances, geometry, and
covariance residuals.  Check exact restart identity and rejection after a
physical atom displacement or operation-set change.

**Step 2: Run RED**

Run checkpoint and coefficient-RT runners.  Expected: FAIL because the exact
fragment symmetry fingerprint is not stored.

**Step 3: Implement V3 fingerprint binding**

Extend existing V3 metadata fields/fingerprints without creating V4.  Keep
one-shot and restart-split observable files byte-identical.

**Step 4: Focused verification and reviews**

Run checkpoint and RT on 1/2/4/8 ranks, route/removal contracts, and
`git diff --check`.  Perform both reviews, resolve all Critical and Important
findings, rerun affected checks, and commit:

```bash
git add src/gs/dc/dg_overlapping_wannier_checkpoint.f90 \
  src/gs/dc/dg_overlapping_wannier_types.f90 \
  src/rt/dg/rt_dg_overlapping_wannier.f90 \
  tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90 \
  tests/dg/run_rt_dg_overlapping_wannier_mpi.py
git commit -m "feat(dg): bind fragment symmetry to V3"
```

### Task 7: Validate ideal and displaced material physics

**Files:**
- Modify: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Add: `tests/dg/data/si64_overlapping_wannier_rt/atom_displaced.dat`
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`
- Modify: `docs/plans/2026-07-31-obsolete-dg-route-removal-results.md`

**Step 1: Add displaced-structure RED**

Require ideal Si64 to recover its exact fragment symmetry, Cartesian response
equivalence, and inversion-compatible odd-harmonic dominance.  Require a fixed
displaced snapshot to report the reduced group or C1 and preserve anisotropy
and allowed even harmonics; do not judge those broken-symmetry signals as
failures.  Both structures must pass field-off stability, `J` versus `dP/dt`,
amplitude scaling diagnostics, and time-step convergence.  Spectra remain
strictly polarization-derived.

**Step 2: Run clean-first parent-prerequisite overlay build**

Create a fresh `git archive HEAD` tree, apply the documented parent-prerequisite
overlay and current task diff, configure MPI, ScaLAPACK, EigenExa, and spglib,
then run:

```bash
cmake --build <overlay-build> --clean-first -j1
```

Expected: SALMON links successfully.  Record binary SHA-256 and build options.

**Step 3: Run the complete focused suite**

Run every overlapping-Wannier metadata, metric, construction, fragment-group,
point-group, operator, physical-matrix, solver, SCF, checkpoint, and RT runner
on 1/2/4/8 ranks.  Run the 32-point-group survey, route/removal contracts,
normal DC LCFO plus EigenExa regression, and `git diff --check`.

**Step 4: Run genuine production matrices**

From fresh accepted GS checkpoints run ideal and displaced Si64 field-off,
Cartesian impulse, weak laser, HHG laser, amplitude comparison, and time-step
halving.  Record manifests, hashes, runtime, memory, symmetry residuals, and
polarization-derived spectra.  Reject forbidden DG-route markers.

**Step 5: Final reviews**

Perform full-branch specification and code-quality reviews.  Resolve every
Critical and Important finding, repeat all affected focused and production
verification, and do not weaken physical gates.

**Step 6: Record and commit**

```bash
git add tests/dg/check_si64_overlapping_wannier_response_hhg.py \
  tests/dg/run_si64_overlapping_wannier_response_hhg.py \
  tests/dg/data/si64_overlapping_wannier_rt/atom_displaced.dat \
  docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md \
  docs/plans/2026-07-31-obsolete-dg-route-removal-results.md
git commit -m "test(dg): validate exact-symmetry response physics"
```

### Task 8: Final branch verification and publication

**Files:**
- Modify: `docs/plans/2026-07-31-obsolete-dg-route-inventory.md`
- Modify: `docs/plans/2026-07-31-obsolete-dg-route-removal-results.md`

**Step 1: Re-run clean-first overlay and all acceptance gates**

Repeat Task 7 Steps 2 through 4 from a new temporary directory.  Compare
hashes and numerical results with the recorded evidence.

**Step 2: Audit retained routes**

Confirm that the only DG route is buffer-supported exact-symmetry
overlapping-Wannier GS to V3 checkpoint to generalized-eigenvalue Exp
coefficient RT.  Confirm normal SALMON and normal DC LCFO plus EigenExa remain
unchanged.

**Step 3: Final reviews and commit**

Repeat full specification and code-quality reviews, resolve all Critical and
Important findings, rerun affected verification, record evidence, and commit:

```bash
git add docs/plans/2026-07-31-obsolete-dg-route-inventory.md \
  docs/plans/2026-07-31-obsolete-dg-route-removal-results.md
git commit -m "test(dg): verify exact fragment symmetry route"
```

**Step 4: Push both remotes**

After verification-before-completion, push the branch to both `origin` and
`upstream` and verify both remote refs equal local HEAD.

