# Direct Pseudo-atomic s+p Complement Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the complete pseudo-atomic s+p excitation shell independent of the fragment empty-state cutoff.

**Architecture:** Append raw buffer-periodic pseudo-atomic orbitals to the fragment eigenstate candidate block and zero-extend the occupied coefficient matrix. Reuse the existing overlapping-Wannier construction and all symmetry, inclusion, localization, and quality gates.

**Tech Stack:** Fortran 2008, MPI, Python source contracts, CMake, clean parent-prerequisite overlay.

---

### Task 1: Add a candidate-window invariance RED fixture

**Files:**

- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Steps:**

1. Construct occupied plus raw complete-shell candidates with two
   different irrelevant empty-state counts.
2. Require equal retained real-space projectors and identical target
   rank.
3. Require the production adapter to append `manifest_values`, compute
   their periodic gradients, and zero-extend occupied coefficients.
4. Run the construction and source fixtures and observe RED.

### Task 2: Augment the production candidate block

**Files:**

- Modify: `src/gs/main_dft.f90`

**Steps:**

1. Preserve the configured fragment-eigenstate count separately.
2. Build the complete s+p manifest on the buffer-periodic box.
3. Allocate an augmented value/gradient block containing eigenstates
   followed by raw s+p orbitals.
4. Compute raw seed gradients with `periodic_box_gradients`.
5. Zero-extend the occupied coefficient rows for the appended shell.
6. Call the existing construction with the augmented rank and raw
   `projection_seed_values`.
7. Run route, projection, construction, solver/density, SCF, row-storage,
   and `git diff --check` focused verification.

### Task 3: Re-run the fixed Si64 window pair

**Files:**

- Modify: `docs/plans/2026-07-27-si64-overlapping-wannier-results.md`

**Steps:**

1. Rebuild the clean parent-prerequisite overlay.
2. Run fresh c192 and c256 buffer-5 2x2x2 rows on 8 MPI ranks and one
   OpenMP thread.
3. Verify checkpoint reuse and forbidden-route absence.
4. Run the unchanged raw-evidence window checker.
5. Perform specification and quality reviews and resolve every
   Critical/Important finding.
6. Accept Task 9 only if the gate passes; otherwise keep Task 10 blocked.
