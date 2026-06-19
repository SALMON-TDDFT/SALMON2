# Global Wannier Owned-Density Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use executing-plans to implement this plan task-by-task.

**Goal:** Build DG-RT density from global Wannier functions owned by fragment-center positions, including cross terms with neighboring Wannier functions, without reprojecting or reorthogonalizing the Wannier basis inside each fragment.

**Architecture:** RT reads the global Wannier transform and owner metadata, constructs local real-space pieces of every relevant global Wannier on each fragment, then forms partial densities by looping over Wannier functions whose centers belong to the local fragment. Cross terms are evaluated once from the owner side and fragment/rank partial densities are summed at the end.

**Tech Stack:** Fortran 90, existing SALMON DG-RT data structures, MPI reductions, Wannier90 checkpoint-derived `wannier90_global_basis.bin`.

---

### Task 1: Load Global Wannier Data In RT

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`

**Steps:**
1. Add arrays for global Wannier metadata: owner fragment, center, spread, transform matrix.
2. Read the full `wannier90_global_basis.bin`, not only metadata.
3. While streaming each fragment `wavefunctions.bin`, accumulate local pieces
   `W_frag(local_basis, iw) = sum_band C_frag(local_basis, band) U(band, iw)`.
4. Verify the file is absent-safe: old inputs still run without global Wannier mode.

### Task 2: Owned-Wannier Density Kernel

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`

**Steps:**
1. Add an opt-in path, e.g. `SALMON_DG_GLOBAL_WANNIER_DENSITY=1`.
2. Build a local list of `iw` with `owner_frag(iw) == ifrag`.
3. Build a cross list for `jw`, initially all Wannier functions for 2x2x2 validation; later restrict by neighboring owner fragment.
4. Compute density from `rho^W_ij W_i^*(r) W_j(r)` with Hermitian symmetry:
   diagonal once, off-diagonal as `2 Re`.
5. Accumulate fragment-local density and reduce.

### Task 3: Validation

**Files:**
- Use existing `/tmp/dg_wannier_trace_diamond` inputs.

**Steps:**
1. Build with EigenExa/Wannier90 enabled.
2. Run GS reuse with OMP=2/MPI=8.
3. Run short RT with global Wannier density enabled.
4. Check that the previous singular owned projection is not used.
5. Check electron count and Pz response.
