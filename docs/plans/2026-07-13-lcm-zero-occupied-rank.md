# LCM Zero-Occupied-Rank Fix Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Prevent Local Chern Marker calculations from copying a dummy orbital on orbital-parallel ranks with zero occupied orbitals.

**Architecture:** Preserve the nonzero work-array extents required by existing linear algebra. Carry the real local occupied-orbital count explicitly into the copy helpers and use it as their loop bound in both scalar and SOI paths.

**Tech Stack:** Fortran 2008, Python source-contract test, CMake MPI/EigenExa build

---

### Task 1: Add a failing source-contract regression test

**Files:**
- Create: `tests/dg/check_lcm_zero_local_occupied.py`
- Test: `tests/dg/check_lcm_zero_local_occupied.py`

**Step 1: Write the failing test**

Check both LCM source files for a copy-helper call that passes `nocc_local`, a
matching dummy argument, and a loop bounded by that argument. Reject use of
`size(zbuf,4)` or `size(zbuf,5)` as the occupied count.

**Step 2: Run test to verify it fails**

Run: `python3 tests/dg/check_lcm_zero_local_occupied.py`

Expected: FAIL because the existing helpers infer one local orbital from the
dummy allocation when `nocc_local=0`.

### Task 2: Pass the actual occupied count to both copy helpers

**Files:**
- Modify: `src/rt/rt_local_chern_marker.f90:134,431-468`
- Modify: `src/rt/rt_local_chern_marker_soi.f90:128,481-511`

**Step 1: Write minimal implementation**

Add `nocc_local` to each call and helper signature.  Assign the helper's local
loop count from the new argument rather than from `size(zbuf,...)`.

**Step 2: Run regression test to verify it passes**

Run: `python3 tests/dg/check_lcm_zero_local_occupied.py`

Expected: PASS.

### Task 3: Verify integration

**Files:**
- Test: `tests/dg/check_lcm_zero_local_occupied.py`

**Step 1: Rebuild**

Run: `cmake --build build-mpi-eigenexa -j 2`

Expected: `salmon` and `test_preparations` build successfully.

**Step 2: Inspect scope**

Run: `git diff --check && git status --short`

Expected: no whitespace errors; only the planned LCM files, regression test,
and plan documents are attributable to this fix, while pre-existing unrelated
changes remain untouched.
