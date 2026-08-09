# Wannier90 MLWF Global-LCFO Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the custom OW occupation-block localization acceptance path with a matrix-only Wannier90 MLWF transform while retaining distributed full-system symmetry validation, buffer-aware fragment ownership, and bounded-memory V3 publication.

**Architecture:** SALMON constructs reciprocal-neighbour and projection matrices from distributed Gamma-real LCFO tiles, gathers only memory-gated band-space matrices to the serial Wannier90 3.1 library coordinator, and scatters the returned transform for distributed application.  SALMON—not Wannier90—canonicalizes the gauge, validates all full-system affine and physical receipts, assigns MLWFs to fragment cores, and publishes V3.

**Tech Stack:** Fortran 2008, MPI, Wannier90 3.1 library mode, ScaLAPACK, EigenExa, BLAS/LAPACK, spglib, Python contract tests, CMake clean overlays.

---

## Mandatory execution discipline

For every task below:

1. use `@test-driven-development` and preserve the genuine RED output before production edits;
2. run the focused tests on MPI 1/2/4/8 where applicable;
3. perform separate specification and code-quality reviews;
4. resolve every Critical and Important finding and rerun affected tests;
5. build from a clean `git archive HEAD` overlay plus only the task diff and explicit parent-prerequisite commits; and
6. commit only reviewed task files, leaving the existing Si64/HHG worktree changes unstaged until the final acceptance task.

This plan begins after memory Task 3B1 has provided cyclic metric/residual and
row-owned symmetry-overlap producers.  Tasks W1--W2 consume those producers.
Task W3 removes the replicated custom-localizer consumer and completes memory
Task 3B2 by adding row-owned post-MLWF group validation.  Only then is memory
Task 3B complete and its Tasks 4--6 may proceed.  This ordering avoids creating
a temporary full-representation adapter solely for a localization path that is
being removed; it does not waive any Task 3B receipt or review.

### Task W1: Add a strict Wannier90 library adapter

**Files:**
- Create: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `cmakefiles/Builder/build_wannier90.cmake`
- Create: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_w90_mpi.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the failing adapter contract**

Add a small Gamma-only two-band fixture that requests an MLWF transform and
asserts finite unitary output, finite centers/spreads, and a non-increasing
Wannier90 gauge-dependent spread.  Add rejection cases for unavailable library
support, complex leakage, malformed dimensions, non-finite input, and a failed
library status.  The route checker must reject any OW V3 path that calls
`localize_dg_occupation_blocks` after global LCFO construction.

**Step 2: Run RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because the adapter and runner do not exist and the production
route still calls the custom localizer.

**Step 3: Implement the minimal adapter**

Expose one SALMON routine that validates dimensions and Gamma reality, writes
only the minimal `.win` control file on the coordinator, calls
`wannier_setup`/`wannier_run`, validates `U`, centers, spreads, and total spread,
and returns an explicit collective status.  Isolate all preprocessor/build
handling in this module.  Build the bundled library with its documented serial
library interface even when SALMON uses MPI; do not invoke `wannier90.x` and do
not serialize wavefunction arrays.

**Step 4: Run focused GREEN verification**

Run the W90 fixture on MPI 1/2/4/8 and the route checker.  Require one library
call, identical canonical scalar results, and collective rejection on every
rank.  Configure a fresh MPI+Wannier90 overlay and run `cmake --build <build>
--clean-first -j4`.

**Step 5: Review and commit**

After specification and quality reviews, resolve all Critical/Important
findings, repeat GREEN and the clean overlay, then commit:

```bash
git add src/gs/dc/dg_overlapping_wannier_w90.f90 src/gs/dc/CMakeLists.txt \
  cmakefiles/Builder/build_wannier90.cmake \
  tests/dg/test_dg_overlapping_wannier_w90_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_w90_mpi.py \
  tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "feat(dg): add Wannier90 MLWF library adapter"
```

### Task W2: Assemble memory-gated MLWF matrices from distributed LCFO tiles

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write dense-reference and memory RED tests**

For a deterministic distributed LCFO basis, compare `M(m,n,b)` and `A(m,w)`
against direct dense Gamma references on MPI 1/2/4/8.  Assert that no rank
gathers `Nglobal_grid*Nwann`, exact overflow-checked coordinator bytes are
reported before allocation, and a deliberately low byte limit fails
collectively without entering Wannier90.

**Step 2: Run RED**

Run the W90 and construction fixtures on MPI 1/2/4/8.  Expected: FAIL because
distributed matrix formation, byte estimates, and limit rejection do not
exist.

**Step 3: Implement tiled matrix assembly**

Use existing physical point IDs and orbital ownership to reduce reciprocal
neighbour overlaps and projection anchors directly into band-space blocks.
Gather only those blocks and eigenvalues to the optimizer rank after the
collective byte gate succeeds.  Account for input matrices, library copies,
output transforms, centers, spreads, and communication buffers with checked
integer arithmetic and update current/peak workspace receipts.

**Step 4: Run focused GREEN verification**

Require dense-reference agreement and identical matrix fingerprints on MPI
1/2/4/8, plus the low-limit rejection and route contract.  Repeat the
MPI+ScaLAPACK+EigenExa+spglib+Wannier90 clean-first overlay build.

**Step 5: Review and commit**

Resolve both reviews and commit only the listed files as:

```bash
git commit -m "feat(dg): assemble bounded-memory MLWF matrices"
```

### Task W3: Apply and canonicalize MLWF transforms under full-system symmetry

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_solver_mpi.f90`
- Modify: `tests/dg/check_dg_fragment_symmetry_production.py`

**Step 1: Write failing physical-gauge tests**

Cover inversion-related centers outside fragment cores, degenerate symmetry
orbits, Gamma-real sign ambiguity, rank-dependent input ownership, and buffer
tails.  Require canonical transforms, occupied projector and density equality,
affine cocycle/metric covariance, center-orbit closure, and one core owner per
MLWF.  Inject symmetry-breaking and ambiguous transforms and require V3-path
rejection.

**Step 2: Run RED**

Run W90, localization, solver, and fragment-symmetry fixtures on MPI 1/2/4/8.
Expected: FAIL because MLWF output is not applied or canonicalized and
`main_dft` still uses the custom localizer.

**Step 3: Implement distributed transform application**

Block-scatter the transform, multiply distributed orbital/grid tiles, order
columns by symmetry orbit and periodic center, fix Gamma-real signs by the
largest canonical core component, and resolve permitted degenerate rotations
against LCFO anchors.  Recompute every SALMON physical receipt after
localization.  Route OW through this implementation and retain the custom
localizer only for isolated legacy/unit-test uses that cannot publish OW V3.

**Step 4: Run focused GREEN verification**

Run W90, construction, symmetry, localization, operator, solver/density, and
SCF fixtures on MPI 1/2/4/8.  Require canonical hashes and physical values to
match across rank counts.  Run the route contracts and a clean-first full
feature overlay.

**Step 5: Review and commit**

Resolve all Critical/Important findings and commit as:

```bash
git commit -m "feat(dg): apply symmetry-validated MLWF transforms"
```

### Task W4: Bind MLWF provenance to V3 and genuine Si64 acceptance

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_gs.in`
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/test_si64_harmonic_morphology.py`
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`

**Step 1: Write publication and physics RED tests**

Require V3 to include Wannier90 version/backend, input-matrix fingerprint,
canonical transform fingerprint, total/invariant/gauge-dependent spreads,
coordinator byte estimate/peak/limit, and post-MLWF symmetry receipts.  Require
ideal undisplaced Si64 with 384 retained/128 occupied states, full-system
inversion, no displacement, and no custom-localizer provenance.

**Step 2: Run RED**

Run checkpoint and route tests.  Expected: FAIL because MLWF provenance is not
serialized or required.  Record the genuine Si64 checker's pre-run RED for its
missing accepted checkpoint/results.

**Step 3: Implement publication gates**

Extend manifest size/digest/broadcast validation and reject absent, zero,
noncanonical, over-limit, or rank-inconsistent MLWF provenance.  Preserve the
generalized-eigenvalue Exp-only V3 reader boundary and the normal DC
LCFO+EigenExa route.

**Step 4: Run complete acceptance**

From committed HEAD plus only reviewed task diffs, configure Release with MPI,
ScaLAPACK, EigenExa, spglib, and Wannier90 enabled; build EigenExa with `-j1`
and SALMON with `--clean-first -j4`.  Run genuine Si64 GS to V3 on eight MPI
ranks/one OpenMP thread.  Only after all memory, MLWF, inversion, covariance,
stationarity, and publication receipts pass, run field-off, impulse LR, and
long-pulse Exp coefficient RT.  Compute spectra from polarization, use current
only as a secondary `dP/dt` check, classify H2/H4 as peak/dip/slope, and produce
the absolute-path semilog HHG figure.

**Step 5: Final reviews, verification, commit, and dual push**

Run every retained fixture on MPI 1/2/4/8, route/removal contracts, genuine
Si64 and morphology checkers, `git diff --check`, and a final committed-HEAD
clean overlay.  Complete final specification and quality reviews and resolve
all Critical/Important findings.  Commit only then, push the branch to both
`origin` and `upstream`, and verify all three commit IDs are identical.
