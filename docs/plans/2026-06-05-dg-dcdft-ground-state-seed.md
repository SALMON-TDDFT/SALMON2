# DG DCDFT Ground-State Seed Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make DGDFT/LCFO ground-state coefficients the primary DG-RT initial state and demote RT-side Flux-SCF to a fallback/polish path.

**Architecture:** `wavefunctions.bin` already distinguishes dense DC-LCFO coefficients from fragment-local identity seeds. Dense seeds should be trusted as the DGDFT ground-state route: build the flux-aware DG Hamiltonian, diagnose stationarity, and enter RT without running the local Flux-SCF loop. Identity seeds keep the old fallback behavior because they are not Flux-inclusive ground-state coefficients.

**Tech Stack:** Fortran DG-RT modules, existing DG fragment metadata, existing residual diagnostics.

---

### Task 1: Route Dense DC-LCFO Seeds as Primary Initial State

**Files:**
- Modify: `src/rt/main_tddft.f90`

**Steps:**
1. After `prepare_fragment_local_eigen_basis`, branch on `dg_frag%identity_seed_coefficients`.
2. For dense seeds, enable flux trace mixing for one Hamiltonian rebuild so coefficient-dependent face traces are consistent with the loaded coefficients.
3. Skip `prepare_fragment_flux_scf_coefficients`.
4. Log that the dense DGDFT/LCFO ground-state seed is the primary path.
5. Run the existing residual check and density mismatch diagnostic.

### Task 2: Keep Identity Seed as Fallback

**Files:**
- Modify: `src/rt/main_tddft.f90`

**Steps:**
1. Run the current local eigen contraction and Flux-SCF loop only when the seed is identity/local.
2. Change logs from "primary initial-state path" to "fallback polish path".
3. Preserve current fatal residual guard.

### Task 3: Verify

**Commands:**
- `git diff --check -- src/rt/main_tddft.f90`
- `cmake --build build-mpi -j 4`

**Expected:** Diff check passes. Local build may stop at the known OpenMP_C CMake issue before compiling Fortran.
