# Symmetry-Closed Spectral Basins Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make periodic spectral basins accept row-owned symmetry maps and return the finest generator-invariant coarsening of the preliminary watershed partition.

**Architecture:** Keep the current distributed score and six-neighbour watershed. Add a replicated point-level union-find whose initial classes are the watershed labels, stream one row-owned generator at a time into one global vector, and propagate equivalence through generator images until no class merge occurs. Rebuild canonical basin labels and independently verify the induced generator permutations.

**Tech Stack:** Fortran 2008, MPI, Python fixture runners, CMake, EigenExa overlay.

## Global Constraints

- Do not allocate or retain an `Npoint x Ngenerator` global map.
- Accept `generator_maps(nlocal, ngenerator)` paired with `row_ids(nlocal)`.
- Bound convergence by at most `Npoint - 1` successful union operations.
- Reduce every rank-local failure collectively before a shape-dependent MPI branch.
- Preserve user-owned Si64 test/input changes and the untracked W90 fixture outputs.

---

### Task 1: Lock the Row-Owned Generator Contract

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90:1392-1612`

**Interfaces:**
- Consumes: `row_ids(nlocal)` and `generator_maps(nlocal, ngenerator)`.
- Produces: unchanged `build_dg_periodic_spectral_basins` API and MPI-invariant basin fingerprint.

- [ ] **Step 1: Keep the failing distributed-map fixture**

Construct the test map from each local row ID:

```fortran
allocate(spectral_generator_maps(size(spectral_row_ids),1))
do p=1,size(spectral_row_ids)
  spectral_generator_maps(p,1)=merge(int(spectral_row_ids(p)),6-int(spectral_row_ids(p)), &
    mod(int(spectral_row_ids(p)),2)==1)
enddo
```

- [ ] **Step 2: Verify the pre-fix failure**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected before implementation: multi-rank failure at `invalid periodic spectral basin contract`.

- [ ] **Step 3: Enforce the row-owned shape and distributed basin action**

Require `size(generator_maps,1)==nlocal`; range-check only local entries. For every final source basin and generator, reduce the minimum and maximum target basin observed across owned rows. Reject unless they agree and form a permutation.

- [ ] **Step 4: Verify MPI 1/2/4/8**

Run the construction fixture and require one common `SPECTRAL_BASIN_FINGERPRINT`.

- [ ] **Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90
git commit -m "fix(dg): accept row-owned spectral basin maps"
```

### Task 2: Add a Flat-Plateau Symmetry RED

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Interfaces:**
- Consumes: the unchanged basin builder API.
- Produces: a fixture that requires the minimum invariant coarsening rather than a global-ID watershed split.

- [ ] **Step 1: Add the flat periodic fixture**

Use a small periodic line/grid whose score contains a tolerance-flat connected region and whose nontrivial translation maps the preliminary ID-tied watershed across more than one basin. Supply rows cyclically across MPI ranks and assert success, expected coarsened basin count, and a complete basin permutation.

- [ ] **Step 2: Verify RED**

Run the construction fixture. Expected: `spectral basins are not closed under the generators`.

- [ ] **Step 3: Commit the RED only after observing failure**

Record the failing output in the implementation notes; do not weaken the assertion.

### Task 3: Implement the Minimum Generator-Invariant Coarsening

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90:1392-1612`

**Interfaces:**
- Consumes: preliminary `global_labels(npoint)` and row-owned generator maps.
- Produces: canonical coarsened `global_labels`, `basin_labels`, `basin_orbit_map`, fingerprint, and checked workspace receipt.

- [ ] **Step 1: Add checked workspace for union-find streaming**

Preflight and collectively allocate integer arrays for union parent/rank or canonical root, one `global_generator(npoint)`, and any relabel scratch. Include all arrays in the checked `MPI_MAX` workspace receipt.

- [ ] **Step 2: Initialize the equivalence relation**

Union every pair of points sharing a preliminary watershed label. Canonicalize roots using the smallest point ID only as the class name.

- [ ] **Step 3: Stream and propagate each generator**

For generator `g`, place `generator_maps(p,g)` at `global_generator(row_ids(p))`, combine collectively, and verify every entry is assigned exactly once and is a permutation. For every current class, union all generator images of its members. Repeat full generator passes while any union occurs.

- [ ] **Step 4: Bound and synchronize convergence**

Count successful global class merges. Reject collectively if the count exceeds `npoint-1`, if an MPI call fails, or if ranks disagree on whether a pass changed the partition.

- [ ] **Step 5: Canonically relabel and retain the independent closure gate**

Order final classes by their minimum global point ID, rebuild local labels, then run the existing min/max target-basin and permutation checks unchanged as a postcondition.

- [ ] **Step 6: Verify GREEN**

Run the construction fixture on MPI 1/2/4/8. Require the flat plateau case, existing separated basin case, out-of-range RED, and common fingerprints all to pass.

- [ ] **Step 7: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90
git commit -m "fix(dg): close spectral basins under symmetry"
```

### Task 4: Regression and Si64 Production Verification

**Files:**
- Verify only: `src/gs/main_dft.f90`
- Verify only: `/tmp/si64-row-owned-basins-20260815/inputfile`

**Interfaces:**
- Consumes: the completed basin builder.
- Produces: evidence that production passes the former spectral-basin failure without memory regression.

- [ ] **Step 1: Run focused regression**

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
SALMON_GP5_OVERLAY=/tmp/salmon-wpw-full \
  python3 tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
git diff --check
```

- [ ] **Step 2: Rebuild the production executable**

```bash
cmake --build /tmp/salmon-wpw-full -j2
```

- [ ] **Step 3: Run Si64 with the existing resource policy**

Use 8 MPI ranks and one BLAS/OpenMP thread per rank. Monitor RSS and stop only if RSS exceeds 3 GiB/rank or available memory falls below 8 GiB.

- [ ] **Step 4: Record the production receipts**

Require successful `periodic spectral basin` completion and capture basin count, workspace, point-cogroup/affine closure, and the next reached stage. Confirm no `Npoint x Ngenerator` allocation appears.

- [ ] **Step 5: Final verification and commit**

Run `git diff --check`, inspect `git status --short`, and stage only intended source/test files. Do not stage user-owned Si64 files or W90 outputs.

