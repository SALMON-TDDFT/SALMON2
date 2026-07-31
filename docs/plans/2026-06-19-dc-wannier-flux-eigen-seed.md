# DC Wannier Flux Eigen Seed Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create a one-time DC post-processing seed that diagonalizes the Flux Hamiltonian in the global Wannier space and lets RT start from that prepared seed without using full-system matrices during time evolution.

**Architecture:** DC writes a compact `wannier_flux_eigen_seed.bin` next to `wannier90_global_basis.bin`. RT reads this seed, reconstructs DG coefficients from the global Wannier components, and then proceeds with the existing local DG/Wannier propagation.

**Tech Stack:** Fortran, SALMON DC-LCFO flux export, Wannier90 checkpoint data, EigenExa/LAPACK eigensolvers, MPI rank-distributed RT coefficient storage.

---

### Task 1: Define The Seed File Contract

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`

**Steps:**
1. Add a stream binary file named `wannier_flux_eigen_seed.bin`.
2. Store magic/version, `num_bands`, `num_wann`, `nstate_tot`, `nspin`, eigenvalues, and the Wannier-to-eigenstate matrix.
3. Treat this file as a DC preparation artifact only. RT must not create or update it.

### Task 2: Write The Seed During DC Wannier Preparation

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Steps:**
1. After `write_wannier90_global_basis_file`, construct the Flux Hamiltonian in the Wannier space.
2. For the current full Gamma case, use the already diagonal Flux-LCFO eigenbasis and the Wannier transform to form the equivalent eigenstate rotation.
3. Write `wannier_flux_eigen_seed.bin`.
4. Log that the seed is for RT initialization only.

### Task 3: Read And Apply The Seed In RT Initialization

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`

**Steps:**
1. After `global_wannier_coef` is built, check for `data_dcdft/total/wannier_flux_eigen_seed.bin`.
2. If present, zero `dg_frag%coef`.
3. Reconstruct DG coefficients as `global_wannier_coef * seed_wannier_to_eigen` for local rows only.
4. Preserve state-block ownership and `esp`.
5. Print norm diagnostics. Expected norm should remain near the LCFO seed norm, not collapse.

### Task 4: Verify

**Files:**
- Test input: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_w90_formal_impulse_trace5_seedfix`

**Steps:**
1. Build with `cmake --build build-mpi-eigenexa-wannier-lib -j 4`.
2. Run short RT with `OMP_NUM_THREADS=2`, `mpirun -np 8`, and `SALMON_DG_COEF_TRACE=1`.
3. Expected: `norm_base` is order `1e2`, not `1e-27`.
4. Then run `dt=0.02, nt=1000` and inspect Pz.
