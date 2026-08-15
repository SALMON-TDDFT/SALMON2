# Spectral-density Wannier Seeds Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build symmetry-ready Wannier trial channels from occupied and automatically equal-count, energy-resolved unoccupied densities without assuming atomic or bond centres.

**Architecture:** Add small collective primitives for degeneracy-safe spectral windows and streamed density descriptors, then add periodic basin/orbit construction and projected basin operators.  Integrate the resulting deterministic channel catalog with Wannier90 only after standalone MPI invariance and rank-completeness tests pass; retain the current post-Wannier path as a validator during transition.

**Tech Stack:** Fortran 2008, MPI, LAPACK, existing SALMON DG construction modules, Python MPI fixture runners, Wannier90 library interface.

---

### Task 1: Degeneracy-safe equal-count spectral windows

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

**Step 1: Write the failing test**

Add a test for `build_dg_equal_count_spectral_windows` using an occupied block and
eight unoccupied levels, with a threefold degeneracy crossing a nominal window
boundary.  Require nonnegative finite weights, per-state partition sum one,
whole-cluster assignment/taper symmetry, and identical fingerprint on MPI
1/2/4/8.  Add REDs for rank-disagreeing eigenvalues, NaN, unsorted energies, and
an invalid window count.

**Step 2: Run test to verify it fails**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: compilation fails because the new public procedure does not exist.

**Step 3: Write the minimal implementation**

Implement a collective routine that:

- agrees dimensions, tolerance, eigenvalues, and occupations before branching;
- detects occupied/unoccupied states from occupations;
- places approximately equal-cardinality boundaries only between
  tolerance-separated eigenvalue clusters;
- applies a compact smoothstep taper between adjacent windows;
- checks partition of unity and finite output;
- emits a deterministic int64 fingerprint and checked workspace receipt.

**Step 4: Run test to verify it passes**

Run the construction MPI fixture and require PASS for 1/2/4/8 ranks.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_construction_mpi.py
git commit -m "feat(dg): build equal-count spectral windows"
```

### Task 2: Stream occupied and unoccupied density descriptors

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write the failing test**

Create two synthetic models: an ionic model whose occupied and first unoccupied
density maxima are on different sites, and a covalent model with a shared
bond-region maximum plus a higher interstitial maximum.  Rotate degenerate
state columns by a complex unitary and require unchanged descriptors and
fingerprint.

**Step 2: Verify RED**

Run the focused MPI fixture; expect the missing density-builder symbol.

**Step 3: Implement**

Accumulate `rho_occ`, the complete `rho_unocc_total`, and one
`rho_unocc(:,q)` at a time from row-distributed state values and spectral weights.
Require `rho_unocc_total=sum_q rho_unocc(:,q)` so spectral windows can never
discard part of a shell. Normalize into hole/electron/shared feature
fields, use scaled norm accumulation, check finite magnitudes before squaring,
bind row IDs and spectral-window provenance, and avoid a
`grid x state x window` allocation.

**Step 4: Verify GREEN**

Run MPI 1/2/4/8 and compare descriptor fingerprints across decompositions.

**Step 5: Commit**

```bash
git commit -am "feat(dg): stream spectral density descriptors"
```

### Task 3: Periodic basin extraction and symmetry closure

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write REDs**

Test periodic maxima at a cell boundary, deterministic watershed ownership,
ionic separated basins, covalent shared basins, and complete orbit generation
under Z4 and Z2xZ2 actions.  Corrupt one symmetry map and require collective
rejection before indexed access. Require basin scores to consume the complete
unoccupied density and all window components together; forbid per-window rank
selection.

**Step 2: Verify RED**

Run the focused fixture and observe the missing basin API failure.

**Step 3: Implement**

Use global physical row IDs and prepared generator maps.  Stream neighbour
features, choose maxima and watershed ties lexicographically by physical ID,
and generate complete basin orbits from group generators.  Return row-owned
basin labels, orbit metadata, a fingerprint, and checked peak bytes.

**Step 4: Verify GREEN**

Run MPI 1/2/4/8 and require identical basin/orbit receipts.

**Step 5: Commit**

```bash
git commit -am "feat(dg): close spectral basins under symmetry"
```

### Task 4: Projected basin operators and channel catalog

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_eigenexa.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py`

**Step 1: Write REDs**

Require invariance under degenerate-state unitary rotations, complete retained
rank, whole unresolved local eigenspace blocks, and rejection of incomplete
basin catalogs and non-Hermitian projected operators.

**Step 2: Verify RED**

Run the EigenExa MPI fixture and observe the missing API failure.

**Step 3: Implement**

Stream one `K_b = Psi^H b Psi` at a time, diagonalize its Hermitian form, select
only tolerance-separated eigenspaces, propagate equal ranks across each basin
orbit, and assemble row-owned trial channels plus block target actions.  Use
checked dimensions, collective LAPACK/EigenExa gates, and one-basin workspace.

**Step 4: Verify GREEN**

Run EigenExa tests on MPI 1/2/4/8 and construction regressions.

**Step 5: Commit**

```bash
git commit -am "feat(dg): build spectral basin channel catalog"
```

### Task 5: Wannier90 native symmetry integration

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write REDs**

Require deterministic projections instead of `random`, distinct
`d_matrix_band` and basin-block `d_matrix_wann`, AMN covariance, and failure
when the channel catalog does not span all retained states.  Add a generator
DMN fixture proving full-group closure or reject generator-only DMN if the
Wannier90 minimizer does not preserve it.

**Step 2: Verify RED**

Run W90 MPI and route fixtures and observe the expected contract failures.

**Step 3: Implement**

Pass the trial channels into `wannier_setup`, emit matching deterministic
projection data, write the target block representation to DMN while retaining
the streamed band representation, and aggregate spectral/basin/channel
fingerprints into the post-gauge provenance chain.

**Step 4: Verify GREEN**

Run W90 and construction fixtures on MPI 1/2/4/8, route checks, and Release
build.

**Step 5: Commit**

```bash
git commit -am "feat(dg): drive Wannier90 from spectral basins"
```

### Task 6: Si64 validation and transition cleanup

**Files:**
- Modify only existing tracked Si64 fixture files required by the validated path
- Modify: `docs/plans/2026-08-15-spectral-density-wannier-seeds-design.md`

**Step 1: Run focused validation**

Run Si8 through the former Si64 stopping point, followed by Si64 through
Wannier90 completion.  Record per-rank peak memory, spectral-window ranks,
basin orbit sizes, final spread, point closure, Gamma defect, and provenance.

**Step 2: Compare old validator**

Keep the existing post-Wannier centre canonicalizer read-only and require it to
report no corrective rotation beyond tolerance.

**Step 3: Run final verification**

Run all focused MPI 1/2/4/8 suites, route and obsolete-route checks,
`git diff --check`, and the Release build.

**Step 4: Document and commit**

Document measured acceptance receipts.  Do not add the user-owned modified
Si64 experimental files or root-level W90 fixture outputs unless explicitly
authorized.
