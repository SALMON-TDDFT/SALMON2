# Wannier90 Native Line-Search Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Run the corrected Si64 spectral seed through Wannier90's native history-free parabolic line search before TDDFT verification.

**Architecture:** Replace the production fixed-step keyword with Wannier90's standard trial-step keyword while keeping `num_cg_steps=0`. Preserve all unit, symmetry, and convergence gates and use a fresh Si64 spectral-seed run with replay export.

**Tech Stack:** Fortran, Wannier90 v3.1.0, Python route contract, MPI/OpenMP

---

### Task 1: Lock the generated input contract

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Change the route test to require `trial_step = 2.0d0`, reject `fixed_step`, and retain the `num_cg_steps = 0` requirement.
2. Run `python3 tests/dg/check_dg_overlapping_wannier_route.py` and require failure on the old fixed-step writer.
3. Change only the generated Wannier90 keyword from fixed to trial step.
4. Rerun the route test and require PASS.
5. Run `git diff --check` on the two files.

### Task 2: Build and verify the executable

**Files:**
- Build: `/Users/otobetoshihito/SALMON-dev/verification/20260818-onepass-production-build`

1. Rebuild SALMON and the linked Wannier90 library.
2. Verify the resulting `.win` contract with the focused route test.

### Task 3: Run the Si64 spectral-seed gate

**Files:**
- Source input: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_scf.in`
- Create result directory under `/Users/otobetoshihito/SALMON-dev/verification/`.

1. Create a fresh run input with `dg_ow_w90_initial_projection='spectral'`, `wannier_num_iter=10000`, MPI 8, and OMP 1.
2. Enable replay export into a fresh directory.
3. Run without a wall-time cutoff.
4. Confirm the generated `.win` contains the approved settings and corrected Bohr geometry.
5. Record Wannier90 convergence, spread/gradient history, peak RSS, symmetry gates, and whether Hybrid-SCF begins.
6. Do not enter TDDFT validation unless Wannier90 and the subsequent basis gates pass.

