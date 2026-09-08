# Spatial Representation Convention Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Establish one end-to-end spatial representation convention and prevent order-one noncovariant Hamiltonians from being published through group averaging.

**Architecture:** A minimal complex MPI fixture derives the orbital representation from a real spatial point permutation, rather than injecting a matrix. The same fixture checks basis, Cartesian-gradient, and directly integrated operator covariance before production code is changed. Once the correct convention is established, production applies that convention consistently and rejects gross covariance failure before averaging.

**Tech Stack:** Fortran 2008, MPI, BLAS, existing DG construction/operator/symmetry test runners.

---

### Task 1: Add the end-to-end convention fixture

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

**Step 1: Write the failing test**

Construct a distributed periodic point set, a nontrivial point permutation, and two complex basis functions whose induced representation is unitary, nonsymmetric, and non-real. Derive gradients with the same finite-difference stencil used in production. Call the production representation builder and require:

- reconstruction of the permuted basis;
- the correct Cartesian gradient covariance;
- covariance of a directly integrated Hermitian scalar operator;
- rejection of the three incorrect transpose/conjugation alternatives.

Emit a rank-independent convention receipt for the runner.

**Step 2: Run the test to verify RED**

Run: `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

Expected: FAIL at the new representation-convention assertion while all earlier assertions still pass.

**Step 3: Record the failing residuals**

Record all four candidate convention residuals in the test output so the correction is selected from evidence rather than assumption.

**Step 4: Commit the RED fixture**

```bash
git add tests/dg/test_dg_overlapping_wannier_construction_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_construction_mpi.py
git commit -m "test: expose spatial representation convention"
```

### Task 2: Correct the representation convention

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90`

**Step 1: Implement the minimal correction**

Change only the orientation at the boundary shown by Task 1: either normalize the representation returned by `assemble_dg_distributed_basis_symmetry_overlap`, or adjust its consumers. Keep one documented contract for the matrix `D` and derive the operator action algebraically from that contract.

**Step 2: Run focused tests to verify GREEN**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_fragment_symmetry_mpi.py
```

Expected: PASS on MPI 1, 2, 4, and 8; the convention receipt must be identical across rank counts.

**Step 3: Commit the correction**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  src/gs/dc/dg_overlapping_wannier_symmetry.f90 \
  tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90
git commit -m "fix: unify spatial representation convention"
```

### Task 3: Prevent forced repair of gross covariance loss

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90`

**Step 1: Write failing rejection tests**

Add fixtures with an order-one gradient covariance defect and an order-one Hamiltonian component defect. Require collective rejection before the averaged matrix is returned.

**Step 2: Verify RED**

Run the fragment-symmetry MPI runner and confirm the new fixtures fail because the current implementation publishes the group average.

**Step 3: Add the minimal gate**

Require raw basis/gradient/operator covariance to be within the documented numerical tolerance before group averaging. Continue to report before/after residuals; use averaging only for roundoff-scale cleanup.

**Step 4: Verify GREEN and regressions**

Run the construction and fragment-symmetry MPI runners on 1, 2, 4, and 8 ranks, the route checker, and `git diff --check`.

**Step 5: Commit**

```bash
git add src/gs/main_dft.f90 src/gs/dc/dg_overlapping_wannier_symmetry.f90 \
  tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90
git commit -m "fix: reject noncovariant DG operators before averaging"
```

### Task 4: Rebuild and rerun Si64

**Files:**
- No source changes expected.

**Step 1: Build the production executable**

Use the existing one-pass production overlay configuration with MPI 8 and OMP 1.

**Step 2: Run Si64 to the post-Wannier Hamiltonian gate**

Capture per-rank RSS and the basis, gradient, kinetic, local, nonlocal, and total Hamiltonian covariance diagnostics.

**Step 3: Evaluate the result**

Success requires raw component residuals at numerical scale and no occupied-boundary cluster created by group averaging. If a component remains large, stop at that component and trace its input payload; do not add another repair layer.
