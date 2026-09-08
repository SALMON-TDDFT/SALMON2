# Spectral Basin Operator Performance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace scalar basin-operator accumulation and dense MPI reduction with tiled Hermitian BLAS and triangular communication.

**Architecture:** Reuse the prepared basin point-index context, pack a bounded number of weighted state columns, accumulate the upper Hermitian triangle with `ZHERK`, then pack and Allreduce only that triangle. Reconstruct the dense Hermitian result in the existing reusable operator buffer.

**Tech Stack:** Fortran 2008, MPI, BLAS `ZHERK`, existing MPI 1/2/4/8 construction fixture.

---

### Task 1: Fix the optimized communication contract with a RED

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write the failing test**

Extend the prepared projector call with an optional reduced-element receipt and
require `Nstate*(Nstate+1)/2`, while retaining equality with the standalone
operator.

**Step 2: Run test to verify it fails**

Run `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py`.
Expected: compilation failure because the receipt argument is absent.

### Task 2: Implement tiled ZHERK and triangular Allreduce

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

**Step 1: Implement the minimal optimized kernel**

Allocate a fixed-width complex tile and packed upper triangle collectively.
Pack `sqrt(weight)*state_values`, call `ZHERK('U','N',...)`, pack upper values,
Allreduce the packed buffer, and reconstruct both triangles.

**Step 2: Verify GREEN**

Run the construction MPI fixture on 1/2/4/8 ranks and require numerical and
fingerprint equality.

**Step 3: Run regression checks**

Run the EigenExa MPI fixture and `git diff --check`.

**Step 4: Commit**

Commit only the construction source, focused fixture, and these two plan files.
