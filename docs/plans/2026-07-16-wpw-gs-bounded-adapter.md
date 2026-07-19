# WPW GS Bounded Adapter Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Connect the fragment-root WPW operator to production DG-DC through bounded GS-neutral H/S callbacks and a transactional LDA/Hartree potential map.

**Architecture:** Import existing LCFO WW blocks only after they are assembled, combine them with distributed sparse WP/PP blocks in a common owner-targeted adapter, and run retained-subspace SCF using owned W/P slices. Rebuild only potential-dependent volume components during SCF and publish state atomically.

**Tech Stack:** Fortran 2008, MPI, SALMON LCFO/LDA/Hartree infrastructure, LAPACK bounded Rayleigh-Ritz, Python source contracts, linked MPI fixtures, CMake.

---

### Task G1: Lock the bounded callback contract

**Files:**
- Create: `tests/dg/test_dg_wpw_gs_bounded_apply_mpi.f90`
- Create: `tests/dg/run_dg_wpw_gs_bounded_apply_mpi.py`
- Create: `tests/dg/check_dg_wpw_gs_bounded_adapter.py`

1. Write a two-rank RED fixture with different owned/support W counts and fragment-root P ownership.
2. Require batched WW/WP/PP H/S action to match a small dense oracle while inputs and outputs contain owned rows only.
3. Add RED cases for stale epoch, nonfinite coefficients, missing stable IDs, invalid fingerprint, and incomplete reverse delivery.
4. Add source assertions rejecting full/global coefficient vectors, global H/S allocation, global owner arrays, and RT-type dependencies.
5. Run the fixture and source contract; require RED for the missing common adapter.

### Task G2: Implement the common coefficient schedule and apply

**Files:**
- Create: `src/common/dg_wpw_bounded_operator.f90`
- Create: `src/common/dg_wpw_owner_exchange.f90`
- Modify: `src/common/CMakeLists.txt`
- Modify: `src/rt/dg/rt_dg_wpw_matrix_free_adapter.f90`
- Test: `tests/dg/test_dg_wpw_gs_bounded_apply_mpi.f90`

1. Define the plain context, compact WW blocks, sparse WP/PP blocks, owned/input ID lists, epoch, and fingerprint.
2. Move/generalize owner-targeted exchange into `src/common` and make RT consume it; forbid common-to-RT dependencies.
3. Preflight fixed neighbor counts collectively; post receives before sends and complete all requests before status reduction.
4. Implement batched H/S apply and global Gram using only owned/support slices.
5. Store deterministic ordering, metric/operator conventions, layout hashes, and component provenance required by Tasks 6-7.
6. Build a candidate context and replace the live context only after every schedule/block validation succeeds.
7. Run the G1 fixture, owner fetch/reduce tests, existing RT adapter fixture, and `git diff --check`; require GREEN.

### Task G3: Import LCFO WW blocks after assembly

**Files:**
- Create: `src/gs/dc/dg_wpw_lcfo_operator_adapter.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Create: `tests/dg/test_dg_wpw_lcfo_ww_import.f90`
- Create: `tests/dg/run_dg_wpw_lcfo_ww_import.py`

1. Write RED fixtures for local WW, canonical cross-face WW, and distinct periodic minus/plus faces.
2. Import stable W IDs from `n_basis/index_basis` and keep fixed kinetic, fixed nonlocal, fixed face, and potential-dependent WW blocks separate without a global fragment loop.
3. Move sparse operator publication and callback binding after `calc_hamiltonian_matrix`; keep geometry/trace preparation earlier.
4. Fingerprint WW component provenance and the overlap identity/metric convention.
5. Run WW import, canonical-face, operator, full build, and `git diff --check`; require GREEN.

### Task G4: Add the transactional real-space potential map

**Files:**
- Create: `src/gs/dc/dg_wpw_lcfo_potential_map.f90`
- Create: `src/gs/dc/dg_wpw_lcfo_volume_operator.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Create: `tests/dg/test_dg_wpw_lcfo_potential_map_mpi.f90`
- Create: `tests/dg/run_dg_wpw_lcfo_potential_map_mpi.py`

1. Write a RED two-rank one-map fixture with normalized windows, occupied coefficients, charge, density-only mixing, LDA/Hartree update, and changed potential-dependent WW/WP/PP volume blocks.
2. Require one tuple `(rho_in,rho_raw,rho_candidate,V_candidate,E_candidate,H_candidate)` and reject cross-epoch inputs and separately mixed potential.
3. Add RED rollback cases for nonfinite density/potential, charge failure, incomplete core ownership, and failures after mixing, energy evaluation, or operator rebuild; require mixer history, residual history, iteration state, and publication epoch to remain unchanged.
4. Extract callable WW and WP/PP volume rebuild kernels and use them for both initial and iterative operator builds.
5. Implement bounded fragment-core to global-grid redistribution with duplicate/missing-ID and charge-preservation checks, then route potential back only to required support points.
6. Preserve fixed kinetic, nonlocal, and face components; replace only potential-dependent volume components.
7. Compute `V_candidate`, SALMON total-energy components, and the candidate operator from the same `rho_candidate` epoch, including Hartree/XC double-counting corrections.
8. Publish mixer history, iteration state, density, potential, energy, and the new operator epoch together after collective validation.
9. Run the map fixture, LDA/Hartree contracts, rank-local quadrature fixtures, full build, and `git diff --check`; require GREEN.

### Task G5: Connect retained-subspace production SCF

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Create: `src/gs/dc/dg_wpw_scf_driver.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`

1. Write a RED two-level communicator fixture proving every fragment nonroot enters each potential/Hartree/LDA phase and terminates on a root failure.
2. Extend the RED route contract to require an actual call chain from `main_dft` through `lcfo_flux` to the bounded SCF driver.
3. Split root algebra into stepwise operations and implement fixed algebra/potential/converged/failed command and epoch broadcasts.
4. Initialize deterministic stable-ID vectors and retain exactly `n_occ+dg_wpw_extra_states` states.
5. Invoke algebra only after the complete operator epoch is bound, and invoke density/Hartree/LDA on every fragment rank.
6. Fail collectively on every convergence gate and publish converged coefficients/density/potential only after fixed-point density and energy reproduction.
7. Run the communicator fixture, matrix-free SCF fixture, and production consumer contract until GREEN.

### Task G6: Scaling verification and scoped review

**Files:**
- Create: `tests/dg/test_dg_wpw_gs_bounded_scaling.py`
- Inspect: all G1-G5 files and Task 5 production files

1. Add a RED scaling fixture varying fragment count at fixed local basis/support degree.
2. Assert per-root owned/support metadata and apply workspace remain bounded; reject O(N) owner metadata and global fragment scans.
3. Run every `check_dg_wpw_*.py`, all Task 5 MPI fixtures, the full MPI/EigenExa build, and `git diff --check`.
4. Review G1-G5 findings-first. Add a minimal RED regression before each valid fix and repeat until no Critical or Important finding remains.
5. Do not commit, push, or create a pull request.
