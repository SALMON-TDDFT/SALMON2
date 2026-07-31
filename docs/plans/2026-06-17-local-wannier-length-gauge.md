# Local Wannier Length-Gauge DG-RT Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a scalable DC-to-DG route for 10,000+ atom systems using internal
localized Wannier-like basis functions, optional PW augmentation, length-gauge
propagation, and polarization output.

**Architecture:** Keep external Wannier90 seed files as a validation route only.
The production route constructs localized atomic/Wannier-like functions inside the
DC-LCFO-Flux export, removes overcomplete directions locally, exports block-sparse
`H0`, `S`, and position/dipole blocks, and lets DG-RT propagate with
`H(t)=H0-E(t).R`. PW augmentation is added as a local orthogonal complement to
cover nonlocal/high-energy excitations without global dense matrices.

**Tech Stack:** Fortran 90/2003, SALMON DC/DG modules, MPI fragment/orbital distribution, BLAS matmul paths, CMake.

---

### Task 1: Add Explicit Namelist Controls

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_gs`

**Steps:**
1. Add `yn_dc_lcfo_local_wannier`, `yn_dg_length_gauge`, and local Wannier/PW numeric controls.
2. Broadcast, log, lowercase/validate the new controls.
3. Keep defaults off, so existing DC/DG runs are unchanged.
4. Build with `cmake --build /tmp/salmon_wannier90_mpi_config -j 4`.

### Task 2: Isolate Projection Primitives

**Files:**
- Create: `src/gs/dc/dc_lcfo_local_wannier.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/lcfo_flux.f90`

**Steps:**
1. Move or duplicate the C:sp3 projection map/value helpers into a small DC module.
2. Add analytic position matrix support for local Gaussian/SP3 projections.
3. Keep external Wannier90 `.amn/.mmn` writer using the same helper routines.
4. Build and verify no behavior change with defaults.

### Task 3: Export Internal Local Wannier Blocks

**Files:**
- Modify: `src/gs/dc/dc_lcfo_local_wannier.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`

**Steps:**
1. Build local projection basis per fragment core/buffer.
2. Compute local `S`, remove directions below `lambda_cut`, and form an orthonormal local basis.
3. Export fragment-local/block-neighbor `basis`, `H0`, `S`, and `R_x/R_y/R_z`.
4. Add a guard that rejects global dense matrix or `.mmn` production for production-scale settings.

### Task 4: Load Length-Gauge Operators in DG-RT

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`
- Modify: `src/rt/dg/rt_dg_fragment_hamiltonian.f90`
- Modify: `src/rt/dg/rt_dg_integrator_derivative.f90`

**Steps:**
1. Add block-sparse storage for `R_x/R_y/R_z`.
2. Load DC-exported local Wannier blocks.
3. Assemble/apply `H(t)=H0-E(t).R` without velocity-gauge `A.p`.
4. Keep the existing velocity route available for current DG-LCFO tests.

### Task 5: Add Polarization Output

**Files:**
- Modify: `src/rt/dg/rt_dg_observables.f90`
- Modify: `src/rt/dg/rt_dg_iteration.f90`
- Modify: `src/io/write.f90`

**Steps:**
1. Compute `P(t)=-Tr[D(t)R]/Omega` from block-sparse density and position matrices.
2. Print compact stdout columns for polarization when length gauge is active.
3. Write `*_rt.data` consistently with stdout.
4. Verify field-on raw data shows small transverse drift and oscillatory longitudinal polarization.

### Task 6: Add PW Orthogonal Complement

**Files:**
- Modify: `src/gs/dc/dc_lcfo_local_wannier.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`

**Steps:**
1. Generate local PW candidates inside each fragment/buffer.
2. Project out the local Wannier subspace and remove small-overlap directions.
3. Export the mixed `W+PW` basis and all `H/S/R` blocks.
4. Validate convergence of optical response versus PW cutoff.
