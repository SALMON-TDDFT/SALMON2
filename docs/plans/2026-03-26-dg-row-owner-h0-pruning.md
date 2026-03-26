# DG Row-Owner H0 Pruning Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Restrict DG RT `H0` block application to row-owner-local Hamiltonian blocks so each rank computes only the rows it ultimately owns.

**Architecture:** Keep the current row-distributed coefficient ownership model and full `nstate_tot` state batch. Add a lightweight per-rank list of owned Hamiltonian block indices, rebuild that list whenever block topology changes, and route `apply_matrix_blocks_batch` through the owned-block list instead of the global block array.

**Tech Stack:** Fortran 90, SALMON DG RT modules, existing block Hamiltonian representation.

---

### Task 1: Add owned-block metadata to DG fragment state

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90`

**Step 1: Add fields for row-owner-local H blocks**

Add allocatable metadata to `type(s_dg_fragment_rt)` for:
- owned block ids
- owned fragment rows

Keep the representation simple:
- `integer, allocatable :: H_local_block_ids(:)`
- `integer, allocatable :: H_local_rows(:)`

**Step 2: Verify no duplicate/conflicting existing fields**

Run:

```bash
rg -n "H_local_block_ids|H_local_rows" /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg
```

Expected: no matches before implementation.

**Step 3: Commit**

```bash
git add /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90
git commit -m "refactor: add DG local H block metadata"
```

### Task 2: Build the row-owner-local H block list

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_hamiltonian.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_hmat_reconstruct.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90`

**Step 1: Write a small helper that derives owned row fragments**

In `rt_dg_fragment_ops.f90`, add a helper that scans fragment rows and records rows whose mapped basis owner is the current rank.

Rule:
- A fragment row is local if its active basis exists and its first active global basis row is owned by `dg_frag%id`.

**Step 2: Build owned block ids from row fragments**

Add a helper that:
- loops over `H_mat_blocks`
- keeps only blocks whose `ifrag_row` is in the owned row set
- stores block indices in `dg_frag%H_local_block_ids`

**Step 3: Call the helper after block creation/reconstruction**

Invoke the helper after:
- initial Hamiltonian block setup
- Hamiltonian reconstruction
- basis-update block rebuild paths

**Step 4: Add minimal diagnostics**

Add a root-only or rank-local debug line guarded by existing style that can report:
- total H blocks
- local H blocks

This should stay short and low-noise.

**Step 5: Commit**

```bash
git add /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90 /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_hamiltonian.f90 /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_hmat_reconstruct.f90 /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90
git commit -m "refactor: build row-owner local DG H block lists"
```

### Task 3: Route H0 block apply through the local block list

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90`

**Step 1: Add a block-list-aware apply routine**

Implement a variant of `apply_matrix_blocks_batch` that accepts:
- `blocks(:)`
- `block_ids(:)`

and loops only over `block_ids`.

Keep semantics unchanged:
- full state batch
- same accumulation into `y`

**Step 2: Switch RT `H0` apply to the local block list**

In `calculate_time_derivative`, replace the global block apply call with the local-list variant.

Behavior must remain:
- `coef_all(1:n_frag,:)` input
- `dcoef_dt_h0(1:n_frag,:)` output
- later write-back still filtered by `coef_owner`

**Step 3: Preserve fallback behavior**

If the local block list is absent, fall back to the full block list so the code remains safe during initialization/debugging.

**Step 4: Commit**

```bash
git add /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90 /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90
git commit -m "perf: restrict DG H0 apply to owner-local blocks"
```

### Task 4: Verify behavior and build

**Files:**
- Test: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build`

**Step 1: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds.

**Step 2: Sanity-check call sites**

Run:

```bash
rg -n "apply_matrix_blocks_batch|H_local_block_ids" /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg
```

Expected:
- RT derivative uses the local-list path
- metadata is initialized in the intended rebuild paths

**Step 3: Review residual risk**

Check for:
- stale local block lists after basis updates
- empty local block lists on valid ranks
- accidental change in ownership/write-back semantics

**Step 4: Commit**

```bash
git add /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-26-dg-row-owner-h0-pruning.md
git commit -m "docs: add DG row-owner H0 pruning plan"
```
