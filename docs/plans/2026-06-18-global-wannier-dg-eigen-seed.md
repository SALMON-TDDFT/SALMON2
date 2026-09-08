# Global Wannier DG Eigen Seed Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Use whole-system Wannier functions as a basis, diagonalize the DG Hamiltonian projected into that basis, and seed RT with the resulting Flux-DG eigenstates.

**Architecture:** DC export writes fragment-local BPW seed files from the whole-system Wannier90 transform. Each fragment owns a balanced subset of Wannier functions, projects local DG Hamiltonian and position matrices into that subset, diagonalizes the local projected Hamiltonian, and stores eigenstate coefficients `C = W U` in the existing BPW seed format so RT starts from DG-Hamiltonian eigenstates rather than raw Wannier functions.

**Tech Stack:** Fortran 90, SALMON DC-LCFO/DG-RT, EigenExa-enabled build, Wannier90 checkpoint and `*_r.dat`.

---

### Task 1: Add a BPW Seed Regression Probe

**Files:**
- Create: `tmp_input_rt_globalbpw_20`
- Test with: `/tmp/dg_w90_fullaa_seedcheck/inputfile_rt_globalbpw_20`

**Steps:**
1. Build SALMON.
2. Run the ecut=1.0 GS reuse/import case.
3. Run the short RT case with `SALMON_DG_BPW_POSITION_MODE=AA_R`.
4. Expected before the fix: RT stops with invalid occupied seed S norm.
5. Expected after the fix: RT passes initialization and writes non-header RT rows.

### Task 2: Diagonalize Projected DG Hamiltonian in BPW Export

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Steps:**
1. In `write_wannier90_global_bpw_seed`, after forming `wcoef` and `h_wann`, diagonalize `h_wann`.
2. Replace stored `wcoef` by `wcoef * eigvec`.
3. Rotate `r_wann`, `v_wann`, `aa_wann`, centers, spreads, and local-norm diagnostics consistently.
4. Store `h_wann` as diagonal eigenvalues in the new basis.

### Task 3: Verify

**Commands:**
- `cmake --build build-mpi-eigenexa-wannier-lib --target salmon -j4`
- `mpirun -np 8 ... < inputfile_gs_cluster_global_ecut1`
- `mpirun -np 8 ... < inputfile_rt_globalbpw_20`

**Success Criteria:**
- Build succeeds.
- BPW seed logs show 64 states per fragment.
- Short RT reaches propagation and produces data rows.
- If RT still fails, the failure must move past initial occupied S-norm validation.
