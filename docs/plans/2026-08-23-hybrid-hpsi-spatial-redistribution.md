# Hybrid Hpsi Spatial Redistribution Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Apply total-system `hpsi` to row-owned overlapping-Wannier tiles when their MPI grid ownership differs from `dc%mg_tot` ownership.

**Architecture:** Add a cached bidirectional ID permutation and `MPI_Alltoallv` value transpose to the existing full-cell module.  Build the schedule once from source and destination physical IDs, reuse it for every 16-orbital tile, and connect the production callback without creating replicated full-cell orbitals.

**Tech Stack:** Fortran 2008, MPI, SALMON real-space `hpsi`, Python route checks.

---

### Task 1: Specify Cross-Layout Redistribution

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py`

**Step 1: Write the failing MPI test**

Add destination IDs owned in contiguous spatial blocks while source IDs use a deliberately permuted cyclic distribution.  Require forward values to match `value(state,id)`, reverse values to restore source order, and a second application to reuse the schedule with different values.

**Step 2: Add collective adverse cases**

Require duplicate source IDs, missing destination IDs, and out-of-range IDs to fail on all ranks without deadlock.

**Step 3: Run the test and verify RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py
```

Expected: compilation fails because the cached redistribution type and routines do not exist.

**Step 4: Commit the failing test**

```bash
git add tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90 tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py
git commit -m "test: specify full-cell cross-layout redistribution"
```

### Task 2: Implement The Cached Bidirectional Schedule

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_full_cell.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90`

**Step 1: Add the public schedule type**

Store communicator/layout fingerprints, forward and reverse counts/displacements,
packing permutations, local extents, initialization state, and workspace receipt.
Do not store orbital values persistently.

**Step 2: Build the schedule collectively**

Validate both global ID sets, determine destination owner and position for each
source ID, exchange metadata once with `MPI_Alltoallv`, and derive the inverse
unpacking permutation.  Use 64-bit extent arithmetic before allocation.

**Step 3: Implement value-only forward and reverse application**

Accept `complex(real64) values(tile_width,local_count)`.  Pack by cached
permutation, perform one `MPI_Alltoallv`, and unpack in caller ordering.  Repeat
symmetrically for the reverse path.

**Step 4: Run MPI 1/2/4/8 and verify GREEN**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py
```

Expected: all rank counts print `PASS`; adverse cases fail collectively inside
the test harness and are accepted as expected failures.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_full_cell.f90 tests/dg/test_dg_overlapping_wannier_full_cell_mpi.f90
git commit -m "feat: cache full-cell spatial redistribution"
```

### Task 3: Connect The Production Hpsi Callback

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the failing route assertions**

Require a persistent redistribution schedule, destination IDs generated in exact
`dc%mg_tot` loop order, forward redistribution before `hpsi`, reverse
redistribution afterward, and removal of the false same-owner rejection.

**Step 2: Run the route checker and verify RED**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: failure because production still writes source-owned values directly
into `dc%mg_tot`.

**Step 3: Implement the minimal production connection**

Build local `mg_tot` physical IDs once.  Initialize or validate the cached
schedule at the callback boundary.  Forward-transpose `tile_in`, fill
`tile_psi` in `mg_tot` loop order, call `hpsi`, extract in the same order, and
reverse-transpose into `tile_out`.  Preserve component modes and nonfinite-value
diagnostics.

**Step 4: Verify focused tests**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py
git diff --check
```

Expected: PASS with no whitespace errors.

**Step 5: Commit**

```bash
git add src/gs/main_dft.f90 tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "fix: redistribute hybrid hpsi tiles by spatial ownership"
```

### Task 4: Verify Si64 Through Hybrid Hamiltonian Assembly

**Files:**
- Runtime artifact only: `verification/20260823-si64-hpsi-spatial-redistribution-mpi8-omp1/`

**Step 1: Build the production executable**

Run:

```bash
cmake --build /Users/otobetoshihito/SALMON-dev/verification/20260818-onepass-production-build -j 4
```

Expected: `Built target salmon`.

**Step 2: Run Si64 with unchanged parallelism**

Use MPI 8 and OMP 1, pipe the existing input file on standard input, and set no
wall-time cutoff.

**Step 3: Check the acceptance gates**

Require no ownership-contract failure, finite `H psi`, Hermitian projected
kinetic/local/nonlocal matrices, and progress beyond the first hybrid
Hamiltonian assembly.  Record wall time and maximum RSS per rank.

**Step 4: Run regression checks**

Run the full-cell MPI suite, route checker, hybrid-SCF route checker, and
`git diff --check` again after the Si64 evidence is collected.

**Step 5: Request code review**

Use `requesting-code-review` to review collective ordering, schedule invalidation,
workspace scaling, and numerical equivalence before declaring the fix complete.
