# All-Electron Local Chern Marker Contract Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make scalar and SOI Local Chern Marker output count all electrons and fail safely for unsupported occupations or ill-conditioned dual overlaps.

**Architecture:** Build an idempotent occupied projector from unweighted orbitals, then apply an explicit electron multiplicity to each marker contribution. Validate the occupation contract before expensive work and propagate overlap-inversion failure collectively across orbital ranks.

**Tech Stack:** Fortran 2008, MPI communication wrappers, Python source-contract tests, CMake MPI/EigenExa build

---

### Task 1: Add failing all-electron occupation-contract tests

**Files:**
- Create: `tests/dg/check_lcm_all_electron_contract.py`
- Test: `tests/dg/check_lcm_all_electron_contract.py`

**Step 1: Write the failing test**

Require explicit scalar multiplicities 2 (`nspin=1`) and 1 (`nspin=2`), SOI
multiplicity 1, occupation validation, use of the multiplicity in the final
marker contribution, and absence of `sqrt(local_occ_w)` pre-Lowdin scaling.

**Step 2: Verify RED**

Run: `python3 tests/dg/check_lcm_all_electron_contract.py`

Expected: FAIL against the current weighted-Lowdin implementation.

### Task 2: Implement the sharp all-electron projector contract

**Files:**
- Modify: `src/rt/rt_local_chern_marker.f90`
- Modify: `src/rt/rt_local_chern_marker_soi.f90`

**Step 1: Add occupation validation and multiplicity selection**

For each occupied list, validate occupations against 2 or 1 as appropriate.
Stop with a diagnostic identifying orbital, k point, spin, actual occupation,
and expected multiplicity when the sharp contract is violated.

**Step 2: Remove pre-Lowdin occupation scaling**

Copy unweighted orbitals and Lowdin-orthonormalize only those orbitals.

**Step 3: Apply electron multiplicity**

Multiply each k/spin marker contribution by the validated multiplicity.

**Step 4: Verify GREEN**

Run: `python3 tests/dg/check_lcm_all_electron_contract.py`

Expected: PASS.

### Task 3: Add failing collective-condition tests

**Files:**
- Create: `tests/dg/check_lcm_collective_condition_guard.py`

**Step 1: Write the failing test**

Require the inverse helper to return success rather than warning-only, require
the rowwise wrapper to broadcast that status on `info%icomm_o`, and require all
orbital ranks to stop before coefficient distribution when `rcond<1e-10`.

**Step 2: Verify RED**

Run: `python3 tests/dg/check_lcm_collective_condition_guard.py`

Expected: FAIL because the current code only prints a warning.

### Task 4: Implement collective fail-closed overlap inversion

**Files:**
- Modify: `src/rt/rt_local_chern_marker.f90`
- Modify: `src/rt/rt_local_chern_marker_soi.f90`

**Step 1: Return inversion status from the root helper**

When `rcond<1e-10`, mark the inversion invalid and do not call `zgetrs`.

**Step 2: Broadcast and stop collectively**

Broadcast the logical status and condition estimate from orbital root on
`info%icomm_o`; make every rank stop before distributing inverse rows.

**Step 3: Verify GREEN**

Run: `python3 tests/dg/check_lcm_collective_condition_guard.py`

Expected: PASS.

### Task 5: Verify the complete LCM change

**Files:**
- Test: `tests/dg/check_lcm_all_electron_contract.py`
- Test: `tests/dg/check_lcm_collective_condition_guard.py`
- Test: `tests/dg/check_lcm_zero_local_occupied.py`

**Step 1: Run focused tests**

Run all three Python checks. Expected: PASS.

**Step 2: Build**

Run: `cmake --build build-mpi-eigenexa -j 2`

Expected: `salmon` and `test_preparations` build successfully.

**Step 3: Inspect scope**

Run: `git diff --check && git status --short`

Expected: no whitespace errors and pre-existing unrelated SAWF changes remain untouched.
