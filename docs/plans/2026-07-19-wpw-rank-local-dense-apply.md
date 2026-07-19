# WPW Rank-Local Dense Apply Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace scalar sparse-entry H/S application with BLAS-backed rank-local bounded dense blocks.

**Architecture:** Sparse WW/WP/PP arrays remain authoritative for checkpoint and provenance. Transactional initialization additionally packs bounded dense caches, and apply uses matrix multiplication after the existing owner exchange.

**Tech Stack:** Fortran 2008, MPI, BLAS-backed intrinsic `matmul`, Python source contracts, MPI dense oracle.

---

### Task 1: Specify dense-cache contracts

**Files:**
- Modify: `tests/dg/check_dg_wpw_bounded_index_cache.py`
- Modify: `tests/dg/test_dg_wpw_gs_bounded_apply_mpi.f90`

1. Assert that H/S dense cache fields exist for WW, WP, and PP.
2. Assert that `apply_blocks` contains block `matmul` operations and no loop over sparse entries.
3. Extend the MPI oracle to exercise multiple vectors and duplicate-coordinate accumulation.
4. Run the tests and confirm RED for missing caches.

### Task 2: Pack caches transactionally

**Files:**
- Modify: `src/common/dg_wpw_bounded_operator.f90`

1. Add six allocatable dense cache fields to `s_dg_wpw_bounded_operator`.
2. Allocate them with `stat=` in the candidate object using int64-checked extents.
3. Pack sparse values through the validated endpoint caches using additive assignment.
4. Collectively reject allocation, index, or nonfinite failures before publishing the candidate.
5. Run the focused source contract and MPI oracle; expect PASS.

### Task 3: Replace scalar application

**Files:**
- Modify: `src/common/dg_wpw_bounded_operator.f90`

1. Keep owner-targeted W/P fetches unchanged.
2. Compute W partial as `WW*required_W + WP*owned_P`.
3. Compute owned P as `conjg(transpose(WP))*required_W + PP*required_P`.
4. Keep W-owner reduction and collective finite validation unchanged.
5. Run the focused MPI oracle and production-operator MPI test; expect PASS.

### Task 4: Benchmark and production preflight

**Files:**
- Modify only if required by a failing test: focused test sources
- Create: a new directory below `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_runs/`

1. Build `build-mpi-eigenexa`.
2. Run the multi-vector focused benchmark.
3. Run Si64 with 8 MPI ranks, one OpenMP thread, and `dg_wpw_scf_max_iter=1` in a new directory.
4. Require all ranks to pass the first H/S applications and reach a finite algebra diagnostic or intentional max-iteration failure.
5. Run focused tests, full build, and `git diff --check`.

No commit, push, or PR step is permitted for this worktree.
