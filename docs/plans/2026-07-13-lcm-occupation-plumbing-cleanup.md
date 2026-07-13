# LCM Occupation Plumbing Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove obsolete local occupation-weight plumbing from both LCM implementations without changing behavior.

**Architecture:** Retain global occupation values only for sharp-contract validation. Simplify the occupied-distribution cache to ownership and orbital-index data, and delete helpers whose sole consumers were removed.

**Tech Stack:** Fortran 2008, Python source-contract tests, CMake MPI/EigenExa build, Git

---

### Task 1: Add the failing cleanup contract

**Files:**
- Create: `tests/dg/check_lcm_occupation_plumbing_cleanup.py`

Write a check that rejects `local_occ_w`, `local_occ_weight`,
`local_occ_index`, and `local_occ_global_io` in both modules; rejects
occupation-weight arguments in `build_occ_distribution_cache`; and requires
`occ_w` to remain in sharp occupation validation. Run it and verify RED.

### Task 2: Remove obsolete plumbing symmetrically

**Files:**
- Modify: `src/rt/rt_local_chern_marker.f90`
- Modify: `src/rt/rt_local_chern_marker_soi.f90`

Remove dead declarations, cache arguments, allocations, assignments,
deallocations, and helper routines. Preserve `local_occ_position`. Run the new
check and verify GREEN.

### Task 3: Verify, commit, and publish

Run the cleanup check and the three existing focused LCM checks. Build with
`cmake --build build-mpi-eigenexa -j 2`, run `git diff --check`, commit only the
LCM cleanup files and plan documents, then push the current branch to `origin`
and `upstream` without force.
