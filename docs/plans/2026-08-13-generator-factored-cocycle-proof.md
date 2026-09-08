# Generator-Factored Point-Cogroup Proof Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Preserve the exact affine cocycle proof while reducing dense retained-space closure products from all point pairs to generator-complete relations.

**Architecture:** Validate the spatial permutation/cocycle identity for every ordered pair, then validate dense retained-space matrices only for deterministic generators acting on every point element from both sides. Prove generator reachability using the existing group-generator routine and publish proof-size receipts.

**Tech Stack:** Fortran 2008, MPI, existing distributed row-owned overlap assembly, Python MPI test runners.

---

### Task 1: Add generator-proof REDs

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write the failing tests**

Add fixtures that require a nontrivial generator subset, distinguish generator matrix relations from non-generator integer map relations, and assert a checked-pair count smaller than the full square.

**Step 2: Run the focused test and verify RED**

Run: `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

Expected: FAIL because the validator does not expose or use a generator-factored proof receipt.

### Task 2: Implement generator-complete streaming proof

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`

**Step 1: Select and collectively validate generators**

Call `select_dg_group_generators(point_product, point_identity, ...)`, build a logical generator mask, and verify full reachability before dense work.

**Step 2: Preserve full integer action validation**

Run map composition for all ordered pairs using the pullback convention and left translation cocycle.

**Step 3: Restrict dense products**

Execute retained-space overlap assembly and matrix multiplication only when the left or right operand is a generator. Count each unique ordered pair exactly once.

**Step 4: Publish receipts**

Return generator count and checked-pair count from the validator and emit them in the Si64 diagnostic line. Keep workspace receipt collective and conservative.

**Step 5: Run focused tests**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
git diff --check
```

Expected: PASS, including MPI 1/2/4/8.

### Task 3: Review and production verification

**Files:**
- Review: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Review: `src/gs/main_dft.f90`
- Review: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Run specification and quality reviews**

Require zero Critical and Important findings.

**Step 2: Build the clean overlay**

Copy only tracked task changes into a fresh overlay and build SALMON.

Expected: build completes successfully.

**Step 3: Run ideal Si64 safely**

Use a persistent result directory, MPI 8 / OMP 1 / BLAS 1, low process priority, and periodic checkpoint/log copying. Stop if thermal instability is observed.

Expected: generator proof reports substantially fewer than 2304 dense pairs and reaches V3 without changing numerical tolerance.

**Step 4: Commit implementation**

Stage only the generator-proof source/tests and commit after all verification passes.

