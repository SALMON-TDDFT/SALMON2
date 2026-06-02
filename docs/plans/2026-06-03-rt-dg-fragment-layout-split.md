# RT-DG Fragment Layout Split Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reduce `rt_dg_fragment.f90` recompilation cost by moving low-risk layout and ownership helpers into a separate module.

**Architecture:** Keep `rt_dg_fragment` as the public facade. Move helper routines that do not need access to private module state into `rt_dg_fragment_layout`, then import and re-export them from `rt_dg_fragment` where current callers expect the old module. This keeps behavior unchanged while allowing independent compilation of layout helpers.

**Tech Stack:** Fortran 90 modules, SALMON CMake source list, local `mpif90 -cpp` single-file compile checks.

---

### Task 1: Create Layout Module

**Files:**
- Create: `src/rt/dg/rt_dg_fragment_layout.f90`
- Modify: `src/rt/CMakeLists.txt`

**Steps:**
1. Create `rt_dg_fragment_layout.f90` with module header and public exports.
2. Move pure layout/rank helper routines:
   - `get_fragment_group_root_rank`
   - `fragment_from_rank_address`
   - `wrap_global_grid_index`
   - `get_fragment_grid_sender_rank`
   - `wrap_fragment_cartesian_index`
   - `cartesian_index_to_fragment`
   - `find_density_grid_owner`
3. Add the new source before `rt_dg_fragment.f90` in `src/rt/CMakeLists.txt`.
4. Compile `rt_dg_fragment_layout.f90` alone with `mpif90 -cpp`.

### Task 2: Rewire rt_dg_fragment

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment.f90`

**Steps:**
1. Add `use rt_dg_fragment_layout, only: ...` for moved routines.
2. Remove the moved routine bodies from `rt_dg_fragment.f90`.
3. Preserve public exports from `rt_dg_fragment` for any external callers that still import them there.
4. Compile `rt_dg_fragment.f90` with and without `USE_SCALAPACK`.

### Task 3: Review Next Split Candidates

**Files:**
- Inspect: `src/rt/dg/rt_dg_fragment.f90`

**Steps:**
1. Measure remaining line count and subroutine list.
2. Identify whether `build_density_grid_owner_maps` and `build_fragment_global_boxes` can move next.
3. Do not move basis preparation in this pass unless Task 1 and 2 are clean.
