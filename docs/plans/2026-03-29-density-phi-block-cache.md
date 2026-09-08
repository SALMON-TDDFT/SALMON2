# Density Phi Block Cache Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the flat `density_phi_cache` path in DG density reconstruction with a block-friendly phi cache that removes the expensive whole-fragment repack and reduces per-block packing overhead.

**Architecture:** Add a density-reconstruction-specific block cache stored on `s_dg_fragment_rt`, build it lazily from `phi_frag`, and consume it directly in `calculate_density_from_fragments`. Keep the existing flat cache for Hamiltonian reconstruction unchanged so the optimization remains scoped and low-risk.

**Tech Stack:** Fortran, OpenMP, existing DG fragment cache/invalidation flow, CMake build.

---

### Task 1: Add block-cache storage to fragment state

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment.f90`
- Modify: `src/rt/dg/rt_dg_fragment_basis_update.f90`

**Step 1: Add new cache fields**
- Add allocatables for `density_phi_block_cache` and metadata/valid flag.

**Step 2: Invalidate/deallocate on basis refresh and fragment teardown**
- Clear the new cache wherever `density_phi_cache` is invalidated or deallocated.

**Step 3: Build**

Run: `make -C build -j2 salmon`
Expected: build succeeds.

### Task 2: Build lazy block cache in density reconstruction

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Introduce cache sizing and lazy build**
- Compute per-fragment block count and allocate `density_phi_block_cache(grid_block_size, nstate_frag, max_blocks, ifrag_count)`.
- Populate it from `phi_frag` once per valid basis snapshot.

**Step 2: Replace flat-cache phi packing**
- Use the block cache instead of `density_phi_cache` inside the block loop.

**Step 3: Add timing for block-cache build**
- Expose the new cost in the density breakdown so the missing project time becomes attributable.

**Step 4: Build**

Run: `make -C build -j2 salmon`
Expected: build succeeds.

### Task 3: Verify no regressions in the touched path

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Sanity-check unchanged numerical path**
- Ensure the block cache zero-pads partial blocks and preserves `phi_blk(1:npt_blk,1:nbf)` semantics.

**Step 2: Rebuild and inspect output**

Run: `make -C build -j2 salmon`
Expected: build succeeds and density timing output now includes the block-cache term.
