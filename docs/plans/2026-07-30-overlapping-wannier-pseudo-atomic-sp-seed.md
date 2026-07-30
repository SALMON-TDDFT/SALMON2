# Pseudo-atomic s+p Seed Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Use pseudo-atomic s+p orbitals, rather than nonlocal projector functions, as the symmetry-complete overlapping-Wannier projection seed.

**Architecture:** Keep the complete-shell manifest and periodic projector evaluator unchanged. Simplify the production adapter so each s/p channel reads `pp%upptbl_ao(:,l,species)` on `pp%rad`; retain `udvtbl` only in the physical pseudopotential operator.

**Tech Stack:** Fortran 2008, MPI, Python source contracts, CMake, clean parent-prerequisite overlay.

---

### Task 1: Fix the production s+p radial source

**Files:**

- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/test_dg_overlapping_wannier_projection_mpi.f90`
- Modify: `src/gs/main_dft.f90`
- Modify if no longer needed: `src/gs/dc/dg_overlapping_wannier_projection.f90`

**Step 1: Write the failing test**

Require every production s/p seed channel to index
`pp%upptbl_ao(:,ll,species)` and reject a branch that supplies
`pp%udvtbl` to `radial_projector`.

**Step 2: Run the RED tests**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_dg_overlapping_wannier_projection_mpi.py
```

Expected: the route contract fails because the adapter still selects
nonlocal projectors.

**Step 3: Implement the minimal adapter change**

For `ll=0,1`, validate `pp%nrps_ao`, `pp%rad`, and
`pp%upptbl_ao(:,ll,species)`, then populate `radial_grid`,
`radial_projector`, and `radial_count` exclusively from those
pseudo-atomic-orbital tables.  Do not alter the Hamiltonian
pseudopotential path.

**Step 4: Run focused verification**

Run the route, projection, construction, solver/density, SCF, and
row-storage tests on their required 1/2/4/8 MPI configurations.  Run
`git diff --check`.

### Task 2: Re-evaluate the Si64 window gate

**Files:**

- Modify: `docs/plans/2026-07-27-si64-overlapping-wannier-results.md`

**Step 1: Build the clean overlay**

Copy only the reviewed feature files onto the clean parent-prerequisite
overlay and rebuild.

**Step 2: Run fresh c192 and c256 rows**

Use 2x2x2 fragments, buffer 5, 48 complete s+p targets, 8 MPI ranks, and
one OpenMP thread.  Run checkpoint reuse for every accepted row.

**Step 3: Run the public checker**

Validate each row and calculate the raw metric and energy window
differences.  Do not weaken or replace the `1e-4` gate.

**Step 4: Review and record**

Perform specification and code-quality reviews, resolve every
Critical/Important finding, and update the result document.  Commit only
if the parent Task 9 gate passes; otherwise keep Task 10 blocked.
