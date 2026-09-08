# Point-Orbit Center Blocks Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace impossible scalar joint diagonalization closure with a deterministic, point-generated system of periodic-center blocks.

**Architecture:** Keep the corrected Jacobi sweep as a seed localizer, then generate each candidate seed's full compact point orbit, cluster orbit images by periodic-position expectation, and accept only a complete orthogonal block decomposition permuted by every point operation. Resolve internal channels with the existing LCFO operator.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, Python test runners.

## Global Constraints

- Keep MPI at 8 ranks for Si64 verification and one BLAS/OpenMP thread per rank.
- Do not relax `dg_ow_symmetry_tolerance` or the final affine proof.
- Do not allocate spatial `N^2` matrices or all-character tensors.
- Preserve unrelated user changes and untracked W90 fixtures.

---

### Task 1: Add point-orbit block RED

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`

**Interfaces:**
- Consumes: `jointly_canonicalize_dg_sector_periodic_position_gauge`
- Produces: a noncommuting point-orbit fixture that requires block closure

- [ ] Add a two-center fixture whose projected position components do not commute but whose supplied point representation exactly swaps the two center subspaces.
- [ ] Assert success, finite stationary objective, point leakage below tolerance, and projector/fingerprint invariance under an input-sector unitary rotation.
- [ ] Run `python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py` and verify the current scalar-frame implementation fails at the point-leakage assertion.
- [ ] Commit the RED independently.

### Task 2: Construct center blocks from a point orbit

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`

**Interfaces:**
- Consumes: stationary `unitary`, `position_tuple`, `point_representations`, `lcfo_operator`
- Produces: reordered `unitary`, `centers`, and a complete block-monomial point action

- [ ] Add checked compact workspace for orbit vectors, orbit centers, cluster labels, block projectors, Gram matrices, and rank metadata.
- [ ] For each deterministic seed column, compute `D_p * seed`, its three expectation phases, and periodic center clusters.
- [ ] Build each cluster span by Hermitian Gram eigendecomposition using the existing tolerance-scaled numerical-rank convention.
- [ ] Accept only mutually orthogonal spans whose ranks sum to `m`; order blocks lexicographically by periodic center.
- [ ] Verify every `D_p` maps each accepted span into exactly one accepted span and reduce the worst leakage collectively.
- [ ] Apply the LCFO eigensolver independently inside each accepted block, preserving exact residual multiplets.
- [ ] Run the focused W90 fixture on MPI 1/2/4/8 and verify GREEN.
- [ ] Commit the implementation and GREEN fixture.

### Task 3: Add adverse contracts and receipts

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

**Interfaces:**
- Consumes: point-orbit block builder from Task 2
- Produces: collective rank-loss/incomplete-orbit/corrupt-action rejection and diagnostic receipts

- [ ] Add REDs for an orbit cluster with deficient span, incomplete total coverage, and a point matrix that leaks between otherwise valid blocks.
- [ ] Add collective failure paths reporting seed, point, block, numerical rank, and leakage without rank-local returns.
- [ ] Include all live compact arrays in the checked workspace receipt and clean partial allocations on collective failure.
- [ ] Run W90 MPI 1/2/4/8, `git diff --check`, and the overlapping-Wannier route checker.
- [ ] Commit adverse-contract coverage.

### Task 4: Production verification

**Files:**
- No production source changes expected beyond Tasks 2-3.

**Interfaces:**
- Consumes: Release `salmon` built from the feature branch
- Produces: Si64 diagnostic evidence through the unchanged final gates

- [ ] Build Release with `cmake --build /Users/otobetoshihito/SALMON-dev/verification/20260814-core-ownership-fix-build --parallel 1`.
- [ ] Run Si64 with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 prterun -n 8` in a fresh verification directory.
- [ ] Record memory per rank, block count/ranks, worst point leakage, final center-orbit residual, and affine-cocycle closure.
- [ ] If Si64 fails, retain the unchanged gate and diagnose the first failing receipt rather than loosening tolerance.
