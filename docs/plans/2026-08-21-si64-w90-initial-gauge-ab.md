# Si64 Wannier90 Initial-Gauge A/B Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a reproducible initial-gauge switch and run the symmetry-identical Si64 random-projection comparison.

**Architecture:** Thread one validated string setting from SALMON input to the Gamma-only Wannier90 setup adapter. Keep spectral initialization as the default and change only `.win` projection generation.

**Tech Stack:** Fortran 2008, MPI, Wannier90 3.1.0, Python route tests.

---

### Task 1: Specify the switch with failing tests

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Add spectral/no-random and random/projections-block assertions.
2. Add an unknown-mode rejection assertion.
3. Run the focused MPI and route tests and confirm RED.

### Task 2: Thread the minimal setting

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Add `dg_ow_w90_initial_projection`, default `spectral`.
2. Validate/broadcast/log only `spectral` or `random`.
3. Pass it to `setup_dg_w90_gamma_library`.
4. Emit the random projections block only for `random`.
5. Run focused MPI 1/2/4/8 and route tests and confirm GREEN.

### Task 3: Run the Si64 comparison

**Files:**
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_scf.in`

1. Set `dg_ow_w90_initial_projection='random'` and retain `wannier_num_iter=200`.
2. Rebuild the SPGLIB/Wannier90/ScaLAPACK/EigenExa executable.
3. Run 8 MPI ranks with `OMP_NUM_THREADS=1` in a fresh result directory.
4. Compare convergence iteration, RMS gradient, spread, symmetry receipts, and peak RSS against r8.

