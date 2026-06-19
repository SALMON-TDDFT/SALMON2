# DG BPW Formal AA_R Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make DG-Wannier length-gauge RT use Wannier90 AA_R only when the global Wannier subspace spans the same BPW target band space used for RT.

**Architecture:** The DC export path writes Wannier90 seed files with a full BPW target band subspace. BPW AA_R projection is accepted only when the global Wannier transform covers all DC-LCFO states used by BPW and has enough Wannier functions to span the local BPW directions. RT then uses `SALMON_DG_BPW_POSITION_MODE=AA_R` as the formal path.

**Tech Stack:** Fortran, SALMON DC-LCFO export, Wannier90 `_r.dat`, DG BPW seed files.

---

### Task 1: Enforce Full BPW Wannier90 Subspace

**Files:**
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`

**Steps:**
1. Keep existing `wannier_num_bands` and `wannier_num_wann` inputs.
2. In DC export, require formal BPW AA_R to satisfy `num_bands_file == dc%nstate_tot`.
3. Require `num_wann_file >= max(nkeep)` at projection time. If not, leave AA_R unavailable and print a clear warning.
4. Recommend inputs with `wannier_num_bands = nstate_tot` and `wannier_num_wann` large enough for BPW.

### Task 2: Preserve AA_R Projection Semantics

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Steps:**
1. Project Wannier90 `_r.dat` through the global transform and local BPW coefficients only when the dimension checks pass.
2. Keep the exported BPW seed format version 3.
3. Add diagnostics for `num_bands_file`, `num_wann_file`, `nkeep`, and `max_abs(AA_R)`.

### Task 3: Validate RT Behavior

**Commands:**
- Build: `cmake --build build-mpi-eigenexa-wannier-lib --target salmon -j4`
- Run short no-field and impulse tests with `SALMON_DG_BPW_POSITION_MODE=AA_R`.

**Expected:**
- Incomplete old seeds fall back to unavailable AA_R or stop when AA_R is explicitly requested.
- Newly generated full BPW seeds load AA_R.
- Transition diagnostics no longer show keV-dominated `fsum_like`.
