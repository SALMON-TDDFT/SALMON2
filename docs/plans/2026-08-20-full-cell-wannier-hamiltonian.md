# Full-Cell Wannier Hamiltonian Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the fragment-partitioned post-Wannier Hamiltonian with a center-owned, spatially distributed, tiled full-cell SALMON Hamiltonian projection.

**Architecture:** Wannier rows remain owned by the ranks containing their centers, while each active orbital tile is distributed over the total-system real-space slabs.  The existing SALMON periodic Hamiltonian applies to one bounded tile at a time, and only projected matrix rows are reduced to their center owners.  The first implementation remains dense and becomes the oracle for a later neighbor-truncated path.

**Tech Stack:** Fortran 2008, SALMON `hpsi`, MPI, EigenExa, Python route checks.

---

### Task 1: Establish the full-cell tiled projection contract

**Files:**
- Create: `src/gs/dc/dg_overlapping_wannier_full_cell.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py`

**Step 1: Write the failing MPI test**

Construct a periodic one-dimensional complex basis embedded in a 3D grid,
distribute spatial rows across MPI ranks, and assign output matrix rows by
orbital center owner.  Apply a known finite-difference kinetic term plus a
nonconstant periodic local potential.  Require tiled output rows to equal an
explicit dense reference for tile widths 1, 2, and 3.

**Step 2: Run the test and verify RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py
```

Expected: FAIL because the full-cell tiled projection primitive is absent.

**Step 3: Implement the minimal primitive**

Add a callback-neutral primitive that accepts distributed basis values,
row-owner IDs, a tile width, and a tile Hamiltonian callback.  It must:

- collectively validate global dimensions and ownership exactly once;
- stream owned orbital tiles in deterministic global-orbital order;
- invoke the callback on distributed spatial values;
- compute local `conjg(basis_i)*Hbasis_j` contractions;
- reduce only the requested output rows to their owners;
- use checked default-integer MPI counts and checked `int64` byte accounting;
- clean all partial allocations after collective failure.

**Step 4: Run GREEN verification**

Run the MPI fixture on 1, 2, 4, and 8 ranks and compare row values and
fingerprints across rank counts and tile widths.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_full_cell.f90 src/gs/dc/CMakeLists.txt \
  tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py
git commit -m "feat: add tiled full-cell Wannier projection"
```

### Task 2: Adapt SALMON's total-system Hamiltonian to bounded Wannier tiles

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_full_cell.f90`
- Modify: `src/gs/main_dft.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90`

**Step 1: Extend the fixture with separate components**

Add kinetic, local, and separable nonlocal projector references.  Require the
primitive to return each component separately and their sum, all Hermitian and
independent of tile width.

**Step 2: Verify RED**

Run the focused MPI test and confirm the component contract is absent.

**Step 3: Add the production adapter**

Create a bounded adapter around SALMON's existing `hpsi` data model.  For each
tile, populate only the total-system orbital fields needed by `hpsi`, use
`dc%mg_tot`, `dc%system_tot`, `dc%info_tot`, `dc%ppg_tot`, the periodic stencil,
and the current total-system local potential, then extract the tile output.
Do not reimplement the stencil or pseudopotential formulas.

If the existing `hpsi` interface cannot accept the current spatial
decomposition directly, add the smallest adapter that converts one tile to its
established `s_orbital` layout.  Do not add a second Hamiltonian implementation.

**Step 4: Verify GREEN**

Run the focused MPI fixture, production overlay build, and `git diff --check`.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_full_cell.f90 src/gs/main_dft.f90 \
  tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90
git commit -m "feat: apply SALMON Hamiltonian to Wannier tiles"
```

### Task 3: Reconstruct the initial full-cell density from occupied Wannier states

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Test: `tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90`

**Step 1: Write RED tests**

Require the initial density to be reconstructed from the canonical Wannier
basis, identity initial coefficients, and retained occupations.  Forbid copying
`dc%rho_tot_s` directly into the one-shot state after the Wannier basis exists.
The MPI fixture must verify electron count and affine covariance.

**Step 2: Verify RED**

Run the route checker and focused fixture; confirm the current raw DC-density
initialization fails.

**Step 3: Implement the minimal reconstruction**

Call the established distributed density reconstruction before the one-shot
Hamiltonian build.  Publish that density to the total-system density/potential
update path and retain the existing electron-count receipt.

**Step 4: Verify GREEN**

Run route, focused MPI 1/2/4/8, construction MPI 1/2/4/8, and overlay build.

**Step 5: Commit**

```bash
git add src/gs/main_dft.f90 tests/dg/check_dg_overlapping_wannier_route.py \
  tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90
git commit -m "fix: initialize full-cell density from Wannier states"
```

### Task 4: Switch production publication to the full-cell operator

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_obsolete_routes.py`

**Step 1: Write the route RED**

Require the production order:

1. final Wannier/Gamma/character gauge;
2. center ownership;
3. occupied density reconstruction;
4. full-cell tiled SALMON Hamiltonian;
5. raw covariance/Hermiticity gates;
6. generalized EigenExa;
7. checkpoint publication.

Forbid `assemble_dg_stitched_weak_operator_rows` and fragment projector
collection in the final production Hamiltonian callback.  Keep those routines
available only for focused legacy tests until later cleanup.

**Step 2: Verify RED**

Run route and obsolete-route checkers and confirm the current stitched call is
detected.

**Step 3: Replace the production callback**

Wire the tiled full-cell rows into `ow_build_hamiltonian`.  Preserve component
receipts, row ownership, operator fingerprints, EigenExa communicator handling,
and checkpoint semantics.  Group averaging may run only after raw component
covariance is already within tolerance.

**Step 4: Verify GREEN**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_dg_overlapping_wannier_obsolete_routes.py
git diff --check
```

Expected: all PASS.

**Step 5: Commit**

```bash
git add src/gs/main_dft.f90 tests/dg/check_dg_overlapping_wannier_route.py \
  tests/dg/check_dg_overlapping_wannier_obsolete_routes.py
git commit -m "fix: publish full-cell Wannier Hamiltonian"
```

### Task 5: Si64 production verification and neighbor-shell evidence

**Files:**
- Modify: `src/gs/main_dft.f90` only if a bounded diagnostic receipt is needed
- Create: `docs/plans/2026-08-20-si64-full-cell-wannier-results.md`

**Step 1: Run Si64 with the agreed resources**

Use exactly 8 MPI ranks and `OMP_NUM_THREADS=1`.  Record per-rank RSS, tile
workspace, runtime, and kinetic/local/nonlocal/total raw covariance.

**Step 2: Verify publication gates**

Require Wannier90 completion, full-cell Hamiltonian assembly, generalized
EigenExa residual, electron count, Hermiticity, and checkpoint publication.

**Step 3: Record neighbor-shell magnitudes**

Without dropping entries, bin `|H_ij|` and `|S_ij|` by periodic center distance
and fragment-neighbor shell.  This is evidence only; do not introduce a cutoff.

**Step 4: Document results**

Record exact commands, diagnostics, peak memory, and the prospective cutoff
envelope in the results document.

**Step 5: Commit**

```bash
git add docs/plans/2026-08-20-si64-full-cell-wannier-results.md
git commit -m "docs: record Si64 full-cell Wannier verification"
```
