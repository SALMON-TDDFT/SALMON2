# Prepared Translation Intertwiner Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build and validate finite-translation element maps once, then reuse them across every character phase construction.

**Architecture:** Add a prepared-action derived type and separate prepare/apply procedures in the construction module. Keep the one-shot procedure as a compatibility wrapper, and update `main_dft` to prepare once outside the character loop.

**Tech Stack:** Fortran 2008, MPI collectives, existing SALMON construction fixture and Python route checker.

---

### Task 1: Add prepared-action API REDs

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Steps:**
1. Import the missing prepared type, prepare routine, and apply routine in the fixture.
2. Add a fixture comparing one-shot and prepared results and checking the preparation receipt.
3. Require production preparation before the character loop and prepared apply inside it.
4. Run both focused tests and confirm they fail because the API/route does not exist.

### Task 2: Implement preparation and prepared apply

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

**Steps:**
1. Add the prepared-action type with allocatable full element maps and provenance receipts.
2. Move character-independent metadata agreement, generator-map gathering, extent checks, allocation, and algebra validation into the prepare routine.
3. Add prepared apply using direct element-map indexing for orbit and covariance loops.
4. Rewrite the one-shot routine as prepare plus apply.
5. Run the construction MPI fixture on 1/2/4/8 and make it pass.
6. Commit the primitive implementation.

### Task 3: Integrate production reuse

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Steps:**
1. Declare one prepared-action object in the overlapping-Wannier adapter.
2. Prepare it after generator maps/catalog construction and before the character loop.
3. Replace per-character one-shot calls with prepared apply.
4. Release the prepared maps after the character schedule.
5. Run the route checker and construction MPI fixture.
6. Commit production integration.

### Task 4: Full verification

**Files:**
- Verify only.

**Steps:**
1. Run `python3 tests/dg/check_dg_overlapping_wannier_route.py`.
2. Run `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py`.
3. Build `/Users/otobetoshihito/SALMON-dev/verification/20260814-core-ownership-fix-build` with `-j 1`.
4. Run `git diff --check` and inspect staged scope.
5. Commit the verified optimization without staging unrelated Si64 test edits or generated W90 files.
