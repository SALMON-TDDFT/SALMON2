# DG Row-Layout Symmetry-Map Reindexing Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Correct the post-Wannier retained-core redistribution and keep all row-indexed symmetry maps consistent with the new physical-ID row layout.

**Architecture:** Add one collective row-layout reindexing primitive beside the existing sparse row exchange. Test it independently with deliberately permuted distributed layouts, then connect it at the single production layout transition: materialize values by physical ID, reindex both maps, and only then publish the new IDs.

**Tech Stack:** Fortran 2008, MPI collectives, SALMON DG construction module, Python route checks, MPI 1/2/4/8 regression runner.

---

### Task 1: Add a failing distributed row-layout regression

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

**Step 1: Import and call the proposed primitive**

Add `reindex_dg_point_maps_between_row_layouts` to the construction-module import list.  Build old and new local physical-ID layouts containing the same global IDs in different per-rank orders.  Build a nontrivial old-layout target permutation from a physical action, call the primitive, and assert that every new-layout source row maps to the new global row carrying the same target physical ID.

**Step 2: Add collective adverse cases**

Call the primitive with a same-rank duplicate/missing physical ID and with an out-of-range old target row.  Require collective rejection.  Keep all rank-dependent fixture branches collective-safe.

**Step 3: Run the focused fixture and verify RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: compilation fails because `reindex_dg_point_maps_between_row_layouts` is not yet exported/implemented.  Record that this is the intended failure.

**Step 4: Commit the RED fixture**

```bash
git add tests/dg/test_dg_overlapping_wannier_construction_mpi.f90 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
git commit -m "test: expose DG symmetry-map row-layout mismatch"
```

### Task 2: Implement the minimal collective reindexing primitive

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Add the public interface**

Export:

```fortran
public::reindex_dg_point_maps_between_row_layouts
```

Use a row-owned interface:

```fortran
subroutine reindex_dg_point_maps_between_row_layouts(comm,old_local_ids,new_local_ids,&
    old_target_rows,new_target_rows,ok,message)
```

where both ID arrays are `integer(int64)`, targets are `integer(int64)`, and the output is allocatable.

**Step 2: Validate metadata before shape-dependent work**

Collectively agree the operation count and global row count.  Guard default-integer and MPI counts with `int64` arithmetic.  Verify old and new physical-ID sets each own every global physical ID exactly once and verify every target row lies in `1:global_row_count`.  Use `STAT=` plus collective allocation consensus.

**Step 3: Convert targets through physical IDs**

Collectively assemble:

```text
old global row -> physical ID
physical ID -> old global row
physical ID -> new global row
old global row -> each target old global row
```

For each new local source physical ID, locate its old source row, obtain the old-layout target row, convert that target to a physical ID, and then to its new-layout global row.  Preserve the operation axis exactly.

**Step 4: Define failure cleanup**

On every collective failure, deallocate partial output and all temporary ownership/lookup arrays.  Return `ok=.false.` with a specific message; never allow one rank to return before peers have reached the same consensus point.

**Step 5: Run the fixture and verify GREEN**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: PASS for MPI 1, 2, 4, and 8, including the reordering and adverse cases.

**Step 6: Commit the primitive**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 tests/dg/test_dg_overlapping_wannier_construction_mpi.f90
git commit -m "fix: reindex DG point maps across row layouts"
```

### Task 3: Connect the corrected transition in production

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py`

**Step 1: Make the route check fail on the old transition**

Require the post-Wannier transition to contain, in order:

1. `materialize_ow_distributed_core_to_buffer` from `global_closed_core` using the old `ow_core_ids` and requested `initial_core_ids`;
2. `reindex_dg_point_maps_between_row_layouts` for `global_symmetry_map`;
3. the same primitive for `fixed_center_symmetry_map`;
4. `ow_core_ids=initial_core_ids`.

Forbid passing `initial_core_ids` to `exchange_dg_point_permuted_orbital_rows`.  Update the memory check to allow this bounded core materialization while still forbidding a full-buffer round trip.

**Step 2: Run route checks and verify RED**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py
```

Expected: the route check fails on the existing erroneous exchange and missing map reindex calls.

**Step 3: Replace the erroneous exchange**

Preserve the old IDs until transition completion.  Materialize `global_closed_core` using physical-ID lookup:

```fortran
call materialize_ow_distributed_core_to_buffer(dc%icomm_tot,global_closed_core,ow_core_ids,&
  initial_core_ids,ow_core_values,ok,message)
```

Reindex both symmetry maps into temporary allocatables.  Move-allocate successful results into the production map variables.  Only then assign the new IDs.  On failure, report the helper message and stop consistently, as this initialization path already treats such failures as fatal.

**Step 4: Run route and focused MPI checks**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
git diff --check
```

Expected: all PASS.

**Step 5: Commit production integration**

```bash
git add src/gs/main_dft.f90 tests/dg/check_dg_overlapping_wannier_route.py tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py
git commit -m "fix: preserve DG symmetry maps across core redistribution"
```

### Task 4: Verify the production build and Si64 root-cause prediction

**Files:**
- Do not modify unrelated Si64 inputs or scripts.
- Inspect runtime logs under the existing verification directory.

**Step 1: Build the production executable**

Use the established production build configuration and confirm the changed sources compile and link.

**Step 2: Run the Si64 case with the established resources**

Keep MPI rank count at 8 and set `OMP_NUM_THREADS=1`.  Do not start a second SALMON job concurrently.  Run through the post-Wannier gradient covariance checkpoint.

**Step 3: Verify the predicted diagnostics**

Confirm:

- post-reorder grid-map/stencil commutator is near numerical tolerance rather than about `1.48`;
- gradient covariance defect is near numerical tolerance rather than about `1.4`;
- values covariance remains small;
- no memory spike or unexpected extra full-grid allocation was introduced.

If the first two do not improve, stop at that checkpoint and retain the logs; do not weaken the physics gates.

**Step 4: Run final verification**

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py
git diff --check
```

Expected: all focused checks PASS and Si64 crosses the previously failing covariance checkpoint.

**Step 5: Commit only any necessary test/diagnostic adjustment**

Do not add generated Wannier output or existing unrelated dirty files.
