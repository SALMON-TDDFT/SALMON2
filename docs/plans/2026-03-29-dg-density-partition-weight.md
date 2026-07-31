# DG Density Partition Weight Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the DG density `1/count` overlap averaging with a cosine-taper partition-of-unity weight so reconstructed charge better respects DC fragment support.

**Architecture:** Build a smooth fragment-local taper weight on each fragment grid point, aggregate its owner-local normalization sum, and use normalized weights on every density accumulation path. Keep the existing owner-distributed communication pattern and only swap the weighting model.

**Tech Stack:** Fortran 90, MPI communication helpers, OpenMP loops, SALMON DG RT data structures

---

### Task 1: Add partition-weight storage to DG fragment state

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment.f90`

**Step 1: Add new fields**

Add fields for:

- fragment-local raw taper weights
- owner-local partition normalization sum
- owner-local normalized partition weight

**Step 2: Remove or deprecate old count-average fields**

Mark old overlap-count weighting fields as replaced:

- `density_weight_local`
- `density_inv_weight_local`

**Step 3: Update deallocation/finalize paths**

Ensure new arrays are freed anywhere density maps are rebuilt or fragment state is finalized.

**Step 4: Build**

Run: `make -C build -j2 salmon`
Expected: successful build

### Task 2: Build cosine taper raw weights in owner-map setup

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment.f90`

**Step 1: Implement 1D cosine taper helper logic**

For each local fragment coordinate:

- return `1` in core region
- return cosine taper in buffer region
- return `0` outside

**Step 2: Form 3D raw fragment weight**

Compute `w = wx * wy * wz` for each fragment-local point while owner maps are built.

**Step 3: Accumulate owner-local normalization sum**

For each point owned by this rank, add raw weights from:

- local fragments
- remote fragments represented in `density_recv_map`

**Step 4: Build normalized owner-local partition weight**

Compute normalized weight as `w / W` using the owner-local sum.

**Step 5: Build**

Run: `make -C build -j2 salmon`
Expected: successful build

### Task 3: Replace density accumulation weighting

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Replace all old weight applications**

Swap old:

- `density_inv_weight_local(...)`

for new normalized partition weights on:

- local accumulation path
- subgroup self path
- recv unpack path

**Step 2: Keep charge budget diagnostics working**

Ensure raw and weighted charge-budget diagnostics still report meaningful totals after the weight swap.

**Step 3: Build**

Run: `make -C build -j2 salmon`
Expected: successful build

### Task 4: Validate charge conservation behavior

**Files:**
- Modify if needed: `src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Run target case**

Run the same DG RT case that currently reports about 3% charge loss.

**Step 2: Inspect logs**

Check:

- `density charge budget`
- `total_charge`
- `rho_scale_factor`
- `norm`

**Step 3: Adjust taper formula only if needed**

If charge drift remains large, tune only the taper/core definition before changing architecture.

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_types.f90 src/rt/dg/rt_dg_fragment.f90 src/rt/dg/rt_dg_density_reconstruct.f90 docs/plans/2026-03-29-dg-density-partition-weight-design.md docs/plans/2026-03-29-dg-density-partition-weight.md
git commit -m "fix: replace dg density count averaging with partition weights"
```
