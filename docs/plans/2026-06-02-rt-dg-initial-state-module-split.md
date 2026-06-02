# RT-DG Initial State Module Split Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reduce rebuild time by moving DG initial-state and distributed eigensolver routines out of the giant `rt_dg_fragment.f90` translation unit.

**Architecture:** Create `rt_dg_initial_state` as a separate Fortran module that depends on `rt_dg_fragment_types` and `rt_dg_fragment_ops`. Keep public procedure names available through `rt_dg_fragment` by use-association, and move numerical helper routines needed by the extracted code into the new module.

**Tech Stack:** Fortran modules, SALMON CMake source list, EigenExa guarded by `USE_EIGENEXA`.

---

### Task 1: Extract Initial-State Routines

- Create `src/rt/dg/rt_dg_initial_state.f90`.
- Move `measure_fragment_initial_surface_residual`, `diagonalize_initial_dg_full_distributed`, `relax_initial_occupied_subspace_block_sparse`, and their helper routines.
- Preserve `#ifdef USE_EIGENEXA` guards.

### Task 2: Wire Module Dependencies

- Add `dg/rt_dg_initial_state.f90` before `dg/rt_dg_fragment.f90` in `src/rt/CMakeLists.txt`.
- Import the extracted public procedures in `rt_dg_fragment.f90` so existing callers can keep using `rt_dg_fragment`.

### Task 3: Verify

- Confirm no duplicate extracted routine definitions remain.
- Run line-length and `git diff --check` checks.
- Build if the local build setup is available.
