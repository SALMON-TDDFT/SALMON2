# Overlapping-Wannier Memory Lifetime Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove unnecessary retained buffer and duplicate transform storage from the Si64 overlapping-Wannier production route.

**Architecture:** First shorten allocatable lifetimes without changing arithmetic. Then independently stream the Gamma transform with one-point work vectors. Static route tests protect memory-shape contracts and focused MPI tests protect numerical behavior.

**Tech Stack:** Fortran 2008, MPI, Python 3 static route checks, EigenExa/Wannier90 focused MPI fixtures.

---

### Task 1: Fix production array lifetimes

**Files:**
- Create: `tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py`
- Modify: `src/gs/main_dft.f90`

**Step 1: Write the failing lifetime checker**

Require the source to deallocate the initial `ow_box_values` and
`ow_box_gradients` after core extraction and before Wannier setup. Require
`move_alloc(global_seed_values,w90_anchors)` instead of a full copy. Require
assembly-only inputs and `lcfo_occupied_core` to be released after last use.

**Step 2: Verify RED**

Run: `python3 tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py`

Expected: FAIL because the buffer arrays and duplicate anchors are retained.

**Step 3: Implement minimal lifetime changes**

Move each deallocation immediately after its proven last read. Replace the
seed-to-anchor assignment with `move_alloc`. Remove the obsolete late
deallocations so every allocatable is released exactly once.

**Step 4: Verify GREEN**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
git diff --check
```

Expected: PASS.

**Step 5: Commit**

Commit only the lifetime checker, `main_dft.f90`, and the two plan documents.

### Task 2: Stream the Gamma transform

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_memory_lifetimes.py`
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

**Step 1: Extend the checker with a failing full-temporary assertion**

Reject allocation of `new_values(nstate,npoint)` and
`new_gradients(3,nstate,npoint)` inside `apply_dg_w90_gamma_transform`.

**Step 2: Verify RED**

Run the lifetime checker and require failure on the existing full temporaries.

**Step 3: Implement point-streamed transformation**

Allocate value and three-axis gradient vectors of length `nstate`. For each
point, compute all transformed components before overwriting that point.
Retain the canonical sign phase and transform/center/spread ordering.

**Step 4: Verify GREEN statically**

Run the lifetime and route checkers plus `git diff --check`.

**Step 5: Verify numerically after the active Si64 job ends**

Run the focused Wannier90 MPI fixture on 1/2/4/8 ranks and compare its direct
frame, gradient, and fingerprint assertions.

**Step 6: Commit**

Commit the streamed transform and checker update separately.

### Task 3: Measure the new Si64 peak

**Files:**
- Reuse: `tests/dg/run_si8_overlapping_wannier_memory.py` monitor helpers
- Create only external evidence under `/Users/otobetoshihito/SALMON-dev/verification/`

**Step 1: Build after the current MPI job exits**

Build the reviewed source with MPI, EigenExa, spglib, and Wannier90 enabled.

**Step 2: Run focused verification**

Run construction/Wannier MPI tests on 1/2/4/8 ranks without another production
job active.

**Step 3: Rerun Si64 MPI8 with the same monitor thresholds**

Compare phase-aligned total and per-rank RSS against
`20260813-si64-observe-mpi8/memory.csv`.

**Step 4: Accept or continue**

Require a bounded peak reduction consistent with removing at least the old
buffer retention. Only then consider more invasive buffer tiling.

