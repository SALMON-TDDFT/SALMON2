# Wannier Symmetry Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove automatic inversion projection from direct AMN imports and prevent unavailable position diagnostics from appearing as zero residuals.

**Architecture:** Keep symmetry enforcement behind the explicit `local_inversion_position` mode. Direct AMN routes retain their gauge-consistent position matrix and run diagnostics only. Separate unavailable position data from a computed residual in log output.

**Tech Stack:** Fortran 2008, Python source-regression tests, CMake, MPI, EigenExa, Wannier90.

---

### Task 1: Add regression expectations

**Files:**
- Modify: `tests/dg/check_direct_amn_global_gauge.py`

1. Require the direct-gauge branch to contain no call to a symmetry projection.
2. Require inversion projection to remain guarded by `local_inversion_position`.
3. Require the global operator diagnostic to print `z_odd_available=F` when the position matrix is unavailable.
4. Run `python3 tests/dg/check_direct_amn_global_gauge.py` and verify that it fails on the current source.

### Task 2: Correct the import control flow

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

1. Replace the direct-gauge automatic projection call with the existing global diagnostic call.
2. Preserve `symmetrize_fragment_wannier_position_import` only in the explicit `local_inversion_position` branch.
3. Run the regression test and verify that the control-flow assertion passes.

### Task 3: Correct unavailable diagnostic output

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

1. Emit the numerical `z_odd_res` only when `position_available` is true.
2. Emit `z_odd_available=F` otherwise, while retaining the Hamiltonian and density permutation diagnostics.
3. Run the regression test and verify that all assertions pass.

### Task 4: Verify build and smoke behavior

**Files:**
- Test: `build-mpi-eigenexa-wannier-lib`
- Test: `build-mpi-eigenexa-wannier-lib/samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac`

1. Run `cmake --build build-mpi-eigenexa-wannier-lib -j 2`.
2. Run the direct-AMN import-only smoke input with 8 MPI ranks and 2 OpenMP threads.
3. Check the log for successful completion, gauge provenance, and truthful position availability output.
4. Review `git diff --check` and the scoped diff before committing.
