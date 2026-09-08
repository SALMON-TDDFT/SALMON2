# DG Continuation Full Review Remediation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove all Important findings from the Task 1-9 full review before zero-field stationarity work.

**Architecture:** Preserve distributed grid and basis ownership from GS through RT. Reuse authoritative periodic-position, DC Hartree/XC, pseudopotential, and energy-decomposition facilities instead of constructing parallel dense substitutes.

**Tech Stack:** Fortran 2008, MPI, SALMON DC/RT infrastructure, standalone Python MPI runners, CMake.

---

Work only in the existing `wpw-s-orthogonal-complement` worktree. Preserve all
unrelated dirty changes and use partial staging for pre-dirty files.

### Task 1: Harden distributed DC density loading

**Files:**
- Modify: `src/gs/dc/dcdft.f90`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Test: the focused divided-DC MPI fixture selected after inspecting existing runners

**Steps:**
1. Add duplicate-ID, missing-ID, and injected count-exchange failure fixtures.
2. Run them and the current contract; require RED for each reviewed defect.
3. Track owner-local multiplicity for every expected point and require exactly one.
4. Check `MPI_Alltoall` and synchronize failure before using counts/displacements; apply the same boundary between both `Alltoallv` calls.
5. Update the source contract to require analyze → closure → freeze and forbid the obsolete builder.
6. Run focused tests at 1/2/4 ranks and commit only Task 1 hunks as `fix(dc): validate distributed continuation density`.

### Task 2: Certify grid catalogs and periodic position provenance

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`
- Modify: `tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`

**Steps:**
1. Add RED fixtures for duplicate/missing/out-of-range grid IDs and a basis crossing the periodic boundary.
2. Compare checkpoint position rows with the existing authoritative periodic DG position builder; require the raw-coordinate implementation to fail.
3. Serialize physical grid extent and the existing position-convention fingerprint.
4. Enforce exactly-once grid ownership and range collectively before publication and after reading.
5. Replace raw coordinate integration with the authoritative wrapped convention without changing legacy readers.
6. Run checkpoint tests at 1/2/4/8 ranks and commit as `fix(dg): certify checkpoint grid and position convention`.

### Task 3: Store meaningful PP and energy provenance

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`

**Steps:**
1. Add RED tests that independently perturb PP identity and each required `s_dft_energy` decomposition field.
2. Define receipts from the existing production PP identity and energy decomposition, not operator fingerprints or arbitrary constants.
3. Hash, serialize, validate, and round-trip those receipts while retaining legacy compatibility.
4. Run checkpoint and GS continuation regressions and commit as `fix(dg): record physical checkpoint provenance`.

### Task 4: Redistribute RT checkpoint shards without dense replication

**Files:**
- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_initialization.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_initialization_mpi.f90`
- Modify: `tests/dg/run_rt_dg_hybrid_initialization_mpi.py`

**Steps:**
1. Add a RED allocation/contract fixture forbidding global dense H/S/component/position arrays during same-rank, expansion, and coalescing reads.
2. Route each serialized row directly to its new cyclic owner by global row ID; rebuild local CSR from stored independent graphs.
3. Route grid shards to physical-grid owners and retain only local basis samples needed for projection.
4. Evaluate residual, orthogonality, Hermiticity, charge, and covariance through distributed sparse actions and scalar reductions.
5. Verify exact payload identity for 1/2/4, 2→4, and 4→2/1; commit as `refactor(rt): redistribute hybrid checkpoint shards`.

### Task 5: Add a genuinely isolated hybrid RT environment

**Files:**
- Modify: `src/rt/initialization_rt.f90`
- Modify: `src/rt/main_tddft.f90`
- Add or modify focused initialization contract tests

**Steps:**
1. Add a RED source/runtime fixture proving the hybrid branch allocates no `spsi_in`, `spsi_out`, or `tpsi` and enters before ordinary restart/orbital setup.
2. Extract the smallest shared grid/PP/Hartree/XC/field initialization into a common helper.
3. Call it from ordinary initialization plus a dedicated hybrid initializer; remove `hybrid_basis_only` and its dummy orbitals.
4. Provide the supported local/semi-local XC path without an ordinary orbital argument.
5. Run ordinary RT build/contracts and commit as `refactor(rt): isolate hybrid continuation initialization`.

### Task 6: Reuse distributed density/Hartree and project owned rows

**Files:**
- Modify: `src/gs/dc/dcdft.f90` or the smallest shared density-potential module
- Modify: `src/rt/main_tddft.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_density_update.f90`
- Add: a production-shaped GS-to-RT MPI smoke runner

**Steps:**
1. Add RED tests forbidding full-grid and `N×N` callback allocations and exercising real grid IDs, PP, Hartree, and XC.
2. Reuse the DC distributed density-to-Hartree entry point with RT-owned grid shards.
3. Project the updated local potential directly into owned operator rows using distributed basis samples; preserve the frozen envelope and structure-keyed exchanges.
4. Require exactly one physical update at time zero and one at each explicit step.
5. Run the production smoke at 1/2/4 ranks plus initialization/checkpoint/length-gauge regressions; commit as `refactor(rt): update hybrid potential distributively`.

### Task 7: Full verification and review

**Steps:**
1. Run every remediation runner named by the 2026-08-30 and 2026-08-27 plans at its supported ranks.
2. Run `cmake --build build-hybrid-commit -j2` and `git diff --check`.
3. Request independent GS/MPI and checkpoint/RT reviews of the complete remediation range.
4. Resolve all Critical/Important findings with new RED tests and task-scoped commits.
5. Stop at the review checkpoint before resuming Task 10.
