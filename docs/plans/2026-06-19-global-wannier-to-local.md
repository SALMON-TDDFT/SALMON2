# Global Wannier To Local Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Convert the global Wannier90 basis into fragment-local Wannier lists without discarding complex phase information.

**Architecture:** Keep `global_wannier_coef` as the source of truth, then build a compact per-local-fragment selected Wannier list based on `global_wannier_owner_frag`. This avoids copying complex global Wannier coefficients into the older real-valued `local_wannier_coef` arrays.

**Tech Stack:** Fortran DG-RT modules, MPI reductions already present in SALMON DG fragment code.

---

### Task 1: Add Local-Selected Global Wannier Storage

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`

**Steps:**
1. Add arrays for selected global Wannier ids, owner fragments, centers, and complex coefficient pieces.
2. Build those arrays after `global_wannier_coef` has been accumulated from `wavefunctions.bin`.
3. Select Wannier functions whose owner fragment is local by default; allow an environment range later.

### Task 2: Use Local-Selected Global Wannier For Observables

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_ops.f90`

**Steps:**
1. Prefer the local-selected global Wannier list for center-polarization when present.
2. Keep the full global projection as a fallback.

### Task 3: Verify

**Commands:**
- `cmake --build build-mpi-eigenexa-wannier-lib -j 4`
- `OMP_NUM_THREADS=2 ... mpirun -np 8 ...`

**Expected:**
- Build succeeds.
- RT prints nonzero `Pz`.
- No NaN/STOP/FATAL.
