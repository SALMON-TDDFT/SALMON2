# WPW Residual-Augmented Matrix-Free Eigensolver Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the production matrix-free WPW algebra converge on a wide-spectrum operator without forming a global dense operator.

**Architecture:** Expand the retained Ritz block with its generalized residual and previous conjugate search block, then solve only the resulting reduced problem. Add a rank-revealing reduced generalized eigensolver so dependent directions are discarded without rejecting a sufficiently large effective subspace.

**Tech Stack:** Fortran 2008, MPI, LAPACK `zheev`, existing WPW H/S and global-Gram callbacks, Python fixture runner.

---

### Task 1: Reproduce the wide-spectrum failure

**Files:**
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Add a second diagonal fixture mode with a wide eigenvalue range and deterministic distributed seeds.
2. Call `run_dg_wpw_matrix_free_algebra_step` and require its lowest retained eigenvalues and residuals.
3. Run `python3 tests/dg/test_dg_wpw_matrix_free_scf.py` and verify RED with algebra `info=40`.

### Task 2: Add the reduced rank-revealing solve

**Files:**
- Modify: `src/common/dg_generalized_algebra.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Add a reduced metric containing a dependent direction to the fixture path.
2. Export a solver that accepts `n`, requested `nev`, H, S, and metric cutoff.
3. Diagonalize S, retain its effective positive subspace, solve the transformed Hermitian H, and return the lowest `nev` pairs.
4. Fail if the effective rank is below `nev`.

### Task 3: Replace Richardson with residual augmentation

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`

1. Initialize a zero previous-search block and form `Z=[Q,R,P]` after the current residual check.
2. Apply H and S to Z through the callbacks and form reduced H/S with `global_gram`.
3. Call the rank-revealing reduced solver for `nretain` states.
4. Set `Q=Z*C`, update `P` from the residual and old-search coefficient blocks, and retain the existing bounded iteration and convergence checks.
5. Run the matrix-free fixture and require GREEN.

### Task 4: Production verification

**Files:**
- Verify: `src/gs/dc/lcfo_flux.f90`
- Verify: `tests/dg/check_dg_wpw_production_consumer.py`

1. Run the focused matrix-free, production-consumer, face-trace, and core-W tests.
2. Run `cmake --build build-mpi-eigenexa -j4`.
3. Run `git diff --check`.
4. Run the two-rank `/tmp/wpw_gs_smoke` case and record the next verified boundary.

No commit step is included because this session explicitly forbids commit/push/PR.
