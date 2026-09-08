# DG Density icomm_frag Project Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Distribute DG density projection across all fragment-subgroup MPI ranks instead of concentrating it on `is_frag_root`.

**Architecture:** Keep the existing packed owner communication on `dg_frag%icomm`, but move the heavy `project` work to all ranks in `dg_frag%icomm_frag`. Root broadcasts the coefficient views needed for one spin/state batch, each subgroup rank computes only its assigned grid blocks, then the subgroup reduces partial density fields before the existing global owner exchange.

**Tech Stack:** Fortran, MPI communicators (`icomm`, `icomm_frag`), OpenMP, SALMON DG RT density reconstruction

---

### Task 1: Add subgroup-project working storage

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`

**Step 1: Add persistent subgroup buffers to the DG fragment type**

Add work arrays for subgroup-broadcast coefficient views and subgroup-local density accumulation only if they are required by the project path.

**Step 2: Initialize and clear the new buffers safely**

Make sure init/deallocation paths in `init_dg_fragment_rt` and final cleanup handle the new arrays without leaking or double-freeing.

**Step 3: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: `Built target salmon`

**Step 4: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_fragment_types.f90 src/rt/dg/rt_dg_fragment.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Add DG subgroup density project buffers"
```

### Task 2: Broadcast raw or mixed coefficient views inside `icomm_frag`

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Add helper for subgroup coefficient broadcast**

Introduce a helper that takes the root-owned coefficient block view and broadcasts it on `dg_frag%icomm_frag`, for both fragment-only and PW/mixed paths.

**Step 2: Keep the helper batch-oriented**

Broadcast only the current spin/state block needed by projection, not the full `nstate_tot` matrix, to cap subgroup memory growth.

**Step 3: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: `Built target salmon`

**Step 4: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_fragment_ops.f90 src/rt/dg/rt_dg_density_reconstruct.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Broadcast DG density batches inside fragment subgroups"
```

### Task 3: Split grid blocks across subgroup ranks

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Replace `is_frag_root`-only project loop with subgroup partitioning**

Assign `igrid0` blocks by `mod(block_index-1, isize_frag) == id_frag`, so every subgroup rank gets a disjoint subset of grid blocks for the same fragment.

**Step 2: Keep owner bookkeeping intact**

Continue using `density_owner_map`, `density_send_slot_map`, and packed owner communication exactly as before, but only for the blocks owned by the current subgroup rank.

**Step 3: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: `Built target salmon`

**Step 4: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_density_reconstruct.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Distribute DG density grid blocks across fragment subgroups"
```

### Task 4: Reduce subgroup partial density before global owner exchange

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Sum subgroup-local `rho` and `rho_s` contributions on `icomm_frag`**

Reduce the local-grid density fields after project completes so root ownership exchange still sees the full fragment-group contribution.

**Step 2: Restrict packed send to fragment roots after subgroup reduction**

Only the fragment root should continue to post `icomm` owner sends, but now from subgroup-reduced density fields.

**Step 3: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: `Built target salmon`

**Step 4: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_density_reconstruct.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Reduce DG density project results inside fragment communicators"
```

### Task 5: Verify the new project path

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Keep only minimal timing or stage logs needed for bring-up**

If temporary logs are needed, keep them fragment-root or max/avg only so they do not mask MPI behavior.

**Step 2: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: `Built target salmon`

**Step 3: Smoke-test a small MPI case**

Run an existing small DG RT case and check that:
- `calculate_density_from_fragments` completes
- total electrons remain finite
- no subgroup rank stalls before the global packed exchange

**Step 4: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -am "Verify DG subgroup density projection path"
```
