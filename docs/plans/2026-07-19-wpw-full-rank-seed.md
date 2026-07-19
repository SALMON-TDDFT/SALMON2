# WPW Full-Rank Deterministic Seed Implementation Plan

> **For Codex:** Use test-driven development and execute this plan task-by-task in the existing protected worktree.

**Goal:** Replace the rank-eight WPW production seed with a deterministic, local, full-rank seed for the 160-state Si64 retained subspace.

**Architecture:** A stateless common initializer maps `(basis namespace, stable ID, state)` to bounded complex coefficients using modular cross terms.  The production LCFO route calls it without changing ownership, communication, metric cutoff, or checkpoint contracts.

**Tech Stack:** Fortran 2008, MPI, LAPACK, Python test runners, CMake.

---

### Task 1: Add the Si64-sized RED regression

**Files:**
- Create: `tests/dg/test_dg_wpw_full_rank_seed.f90`
- Create: `tests/dg/run_dg_wpw_full_rank_seed.py`

1. Construct eight rank-local W/P ID sets with at least 160 global degrees of freedom.
2. Request 160 deterministic retained-state columns from the production initializer.
3. Form the global Euclidean Gram matrix and diagonalize it.
4. Require all 160 eigenvalues to exceed the production relative metric cutoff.
5. Run the test before implementation and confirm it fails because the initializer is absent.

### Task 2: Implement the deterministic seed

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`

1. Add a public initializer with strict shape and stable-ID validation.
2. Generate real and imaginary components from independently salted modular cross terms.
3. Replace the separable trigonometric loops in `run_wpw_production_scf` with the initializer.
4. Run the new regression and confirm all 160 retained directions pass the cutoff.

### Task 3: Focused verification

**Files:**
- Inspect: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Inspect: `tests/dg/run_dg_wpw_gs_bounded_apply_mpi.py`

1. Run the new seed regression.
2. Run matrix-free wide-spectrum and bounded H/S MPI fixtures.
3. Run production source contracts and `git diff --check`.
4. Build `build-mpi-eigenexa` with four jobs.

### Task 4: Fresh Si64 preflight and checkpoint retry

**Files:**
- Modify only the experiment manifest and fresh run-directory provenance files.

1. Create a new run directory; never overwrite or reuse fix4 data.
2. Confirm ranks, threads, memory estimate, free space, input hashes, and pseudo hash.
3. Run a bounded preflight and require routing completion plus at least one successful algebra step.
4. Only then launch a fresh production GS checkpoint calculation.
5. Do not start field-off RT unless GS convergence, finite charge/energy, checkpoint identity, and provenance all pass.

No commit, push, or pull request is allowed for this plan.
