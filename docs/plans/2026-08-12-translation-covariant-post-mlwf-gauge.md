# Translation-Covariant Post-MLWF Gauge Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Restore exact atomic-translation covariance after Wannier90 while retaining a bounded-memory, buffer-supported global Wannier basis for V3 checkpoint and Exp RT.

**Architecture:** Keep the fixed-center DMN and Wannier90 MLWF pass, then discover complete translation orbits and apply the smallest orbit-local polar gauge correction.  Stream one symmetry action at a time and validate the 48 point-cogroup representatives with their translation cocycle instead of retaining 1536 dense matrices.

**Tech Stack:** Fortran 2008, MPI, OpenMP, BLAS/LAPACK, EigenExa, Spglib, Wannier90 library, Python contract runners.

---

### Task 1: Define deterministic translation orbits

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the failing tests**

Add a nontrivial translation fixture with several complete equal-size orbits,
permuted input orbital numbers, periodic centers across a cell boundary, and
MPI 1/2/4/8 ownership.  Require a canonical representative, member-to-group
operation mapping, inverse mapping, and rank-independent fingerprint.  Reject
incomplete orbits, duplicate members, inconsistent translation products,
nonfinite centers, integer overflow, and an orbit split by the fixed retained
rank.

**Step 2: Run RED**

```sh
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: FAIL because the translation-orbit API is absent.

**Step 3: Implement the minimal orbit contract**

Implement `build_dg_translation_orbit_manifest`.  Use the validated translation
product table and periodic center matching.  Use deterministic bipartite
matching and canonical lexicographic representatives; do not depend on MPI rank
or Wannier90 numbering.  Store only `O(Nwann)` metadata.

**Step 4: Focused verification**

Run construction MPI 1/2/4/8, route, obsolete-route, and `git diff --check`.
Require identical manifest fingerprints on every rank.

**Step 5: Reviews and remediation**

Perform a specification review of group-action orientation, periodic matching,
arbitrary fragments, and incomplete-orbit rejection.  Perform a code-quality
review of overflow, collectives, allocation failures, recursion bounds, and
determinism.  Add a new RED for every Critical/Important finding and resolve it.

**Step 6: Clean overlay and commit**

Overlay the committed parent into `/tmp/salmon-oa5b-overlay/source`, run
`cmake --build /tmp/salmon-oa5b-overlay/build --clean-first -j1`, rerun focused
tests, and commit only Task 1 files.

### Task 2: Build the bounded-memory orbit polar correction

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

**Step 1: Write the failing tests**

Create complex-gauge translated orbit blocks whose unconstrained basis is
orthonormal but not translation covariant.  Require the correction to recover
the known gauge, preserve density and rank, satisfy translation products, and
choose the closest unitary polar factor.  Include singular cross overlaps,
split singular-value clusters, nonfinite data, and an intentionally corrupted
translation product.

**Step 2: Run RED**

Run the construction fixture on MPI 1/2/4/8.  Expected: FAIL because the
post-MLWF correction API is absent.

**Step 3: Implement the minimal correction**

Implement row-owned, one-operation-at-a-time cross-overlap assembly and
orbit-local SVD/polar alignment.  Apply a final orbit metric polar correction.
Expose pre/post density defect, orthogonality, translation identity/unitarity/
product residuals, correction norm, spread proxy change, fingerprints, and
peak bytes.  Never allocate `Ntranslation*Nwann^2` or real-space orbit copies.

**Step 4: Focused verification and scaling**

Run MPI 1/2/4/8 and require invariant spectra/fingerprints, decreasing per-rank
large-array bytes, bounded root memory, Gamma-real output when input is real,
and exact adverse-case rejection.

**Step 5: Reviews, clean overlay, and commit**

Review polar orientation, phase convention, degenerate singular subspaces,
metric use, density invariance, MPI count overflow, and error agreement.
Resolve all Critical/Important findings with REDs.  Run the clean-first
full-feature overlay and commit Task 2.

### Task 3: Integrate correction after Wannier90

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_w90_mpi.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write production-route REDs**

Require `Wannier90 -> translation manifest -> orbit polar correction -> full
affine proof -> center redistribution -> V3 checkpoint` in that order.  Require
the correction to act identically on core/buffer values, gradients, transforms,
operator provenance, and fingerprints.  Forbid relaxed center tolerances and
all-1536 dense representation storage.

**Step 2: Run RED**

Run W90 MPI 1/2/4/8 and route checkers.  Expected: FAIL because production
currently checks centers immediately after unconstrained Wannier90.

**Step 3: Implement production integration**

Retain the Wannier90 MLWF transform as the initial gauge.  Construct and apply
the translation correction before recomputing centers.  Then validate the 32
translations and the 48 point representatives with the existing cocycle.
Publish pre/post spread, correction norm, center-orbit residual, density drift,
affine residual, workspace peak, and fingerprint diagnostics.

**Step 4: Focused verification**

Run fragment symmetry, construction, EigenExa, checkpoint, W90, route,
obsolete-route, and Exp-only RT fixtures on MPI 1/2/4/8 plus DMN S3 tests and
`git diff --check`.

**Step 5: Reviews, clean overlay, and commit**

Review physical meaning of the corrected basis, translation versus spatially
varying RT fields, buffer gradients, center ownership with external symmetry
centers, normal DC isolation, and checkpoint provenance.  Resolve all
Critical/Important findings.  Run the full-feature clean-first overlay and
commit Task 3.

### Task 4: Genuine Si64 acceptance through V3

**Files:**
- Modify only if a genuine defect is found in Task 1-3 production files/tests.
- Record evidence in the existing diagnostic log/checkpoint acceptance path;
  do not commit generated run products.

**Step 1: Run the ideal, undisplaced Si64 GS**

Use MPI 8, OMP 1, the successful long-cell acceptance input, and the clean
overlay binary.  Do not add phonons or displaced atoms.

**Step 2: Require the physical and structural gates**

Require rank 384, translation and point-cogroup closure, density/electron-count
invariance, fixed-center DMN publication, Wannier90 completion, bounded spread
increase, periodic center orbits, deterministic fragment redistribution, and
V3 checkpoint publication.  Record peak RSS; do not treat small-cell energies
or harmonic intensities as quantitatively converged physics.

**Step 3: Diagnose any failure with a new RED**

Use systematic debugging.  Do not relax symmetry, rank, density, or center
tolerances.  Add the smallest reproducing fixture, implement the minimal fix,
and repeat focused verification and both reviews.

**Step 4: Final verification and commit**

Run all focused MPI 1/2/4/8 fixtures, DMN/Wannier90 integration, normal DC
LCFO+EigenExa isolation, Exp-only RT, `git diff --check`, and a final
clean-first full-feature build.  Commit any acceptance correction separately.

**Step 5: Continue to response physics only after V3**

After V3 succeeds, run polarization-primary linear response and long-pulse HHG
with Exp propagation.  Plot polarization-derived HHG on a semi-log scale and
classify even orders as peaks or dips.  Current acceptance contains no lattice
vibration, and the small system is qualitative rather than quantitatively
converged.
