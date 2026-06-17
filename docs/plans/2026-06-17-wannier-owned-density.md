# Wannier-Owned DG Density Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a DG-RT density path where Wannier functions are owned by their centers and their full buffered real-space leakage contributes through MPI summation.

**Architecture:** DC export stores local Wannier centers and coefficients. RT reads those fields, can identify the owner fragment for each Wannier function, and reconstructs density by forming real-space Wannier amplitudes from buffered fragment basis functions. Existing fragment-basis density remains available as a comparison path while the Wannier path is validated.

**Tech Stack:** Fortran 90, SALMON DG-RT modules, MPI sparse density exchange, existing `local_wannier_basis.bin` binary format.

---

### Task 1: Export and Read Wannier Centers

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`

**Steps:**
1. Extend `local_wannier_basis.bin` version and append `local_wannier_center(3,nkeep,nspin,ifrag_local)`.
2. Compute centers as the diagonal of `local_wannier_r(:,iw,iw,ispin,i_local)`.
3. Read centers in RT and store them in `s_dg_fragment_rt`.
4. Build and run SALMON.

### Task 2: Owner-By-Center Metadata

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`
- Modify: `src/rt/dg/rt_dg_fragment_layout.f90`

**Steps:**
1. Add `local_wannier_owner_fragment` and `local_wannier_owned` metadata.
2. Map each center to a fragment core using periodic wrapped coordinates.
3. Emit a short diagnostic: owned/total Wannier count per rank.
4. Build and run a GS/RT smoke.

### Task 3: Wannier Density Reconstruction Prototype

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`

**Steps:**
1. Add an environment switch `SALMON_DG_DENSITY_WANNIER=1`.
2. For each owner-local Wannier function, reconstruct buffered real-space values from `local_wannier_coef(:,iw)` and `phi_frag`.
3. Project occupied coefficients to the local Wannier subspace and accumulate `rho += occ*|psi_w(r)|^2`.
4. Reuse existing density send/recv accumulation to sum leakage contributions to grid owners.
5. Print electron count diagnostics for fragment-basis versus Wannier-density paths.

### Task 4: Validate Physics Signals

**Files:**
- Test/run in `/tmp/dg_dc_wannier_flux_blockdiag_omp2`

**Steps:**
1. Re-run DC-GS and confirm block-diag Flux eigenstate export.
2. Re-run RT `nt=10` and confirm `P(0)=0`.
3. Check whether z impulse produces dominant `Pz` rather than comparable `Px/Py`.
4. If `Px/Py` remain large, inspect local Wannier center axes and projection matrix ordering.
