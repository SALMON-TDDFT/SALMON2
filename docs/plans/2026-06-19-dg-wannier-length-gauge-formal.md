# DG-Wannier Length-Gauge Formal Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement the length-gauge DG-Wannier propagation route described in the user specification, while preserving DG neighbor communication and weak scaling.

**Architecture:** Treat the Wannier coefficients `C_{aI}(t)` as the propagation variable.  Store only fragment-local blocks and face-neighbor blocks.  Split the position operator into a local diagonal reference plus a gauge-consistent internal matrix, and keep DG surface correction matrices separate from ordinary Wannier AA_R blocks.

**Tech Stack:** Fortran DG-RT modules, SALMON DC-LCFO export, Wannier90 `_r.dat`, MPI face-neighbor exchange, BLAS/LAPACK or existing EigenExa wrappers where needed.

---

## Review Of The Supplied Specification

The design is acceptable with the following clarifications.

1. `R_I` must be a local reference origin for fragment `I`, not the whole position operator.  The diagonal offset from the actual Wannier center is part of `xi_local`.

2. `xi_local(I)` and `xi_flux(IJ)` must be stored in the same Wannier gauge as the coefficients.  Mixing center-only position data with Wannier90 AA_R is not allowed for production.

3. The off-diagonal neighbor position used in the external-field Hamiltonian must be the DG boundary correction `xi_flux(IJ)`, not an unrestricted full-system AA_R block.  Full AA_R may be used only to construct or validate the local/neighbor blocks, then it must be discarded or cut.

4. `H_flux(IJ)` and `xi_flux(IJ)` are different operators.  `H_flux(IJ)` is the field-free DG surface Hamiltonian block.  `xi_flux(IJ)` is the length-gauge boundary-position correction.  Applying both is not double counting because one belongs to `H0` and the other to `-E(t).r`, but `xi_flux(IJ)` must not include ordinary volume AA_R terms already present in `xi_local`.

5. The static field-free spectrum must be checked before RT.  If the DG flux Hamiltonian collapses the gap, RT cannot be physically meaningful even with a correct length-gauge position matrix.

6. The production route must not keep a full global density matrix, full global position matrix, or all-to-all coefficient exchange during time propagation.

## Required Data Model

Add a new formal storage layer rather than overloading the current BPW/global-Wannier temporary arrays.

### Fragment-Local Blocks

For each local fragment `I`:

- `dg_wannier_n(I)`
- `dg_wannier_ids(:,I)`
- `dg_wannier_owner_frag(:,I)`
- `dg_wannier_center_ref(3,I)`
- `dg_wannier_h0_local(a,b,I)`
- `dg_wannier_xi_local(3,a,b,I)`
- `dg_wannier_coef_basis(mu,a,spin,I)` for density reconstruction and optional projection.
- `dg_wannier_c(a,state,spin,I)` for propagation.

`xi_local` is defined as:

```text
<w_aI | r - R_I | w_bI>
```

so the diagonal center offset is retained in `xi_local`.

### Face-Neighbor Blocks

For each directed face-neighbor pair `(I,J)` owned by the row fragment:

- `dg_wannier_h_flux(I,J;a,b)`
- `dg_wannier_xi_flux(I,J;3,a,b)`
- compact send/receive maps for neighbor `C_{bJ}`.

`xi_flux` is defined as the DG boundary correction:

```text
<w_aI | r - R_IJ | w_bJ>_{boundary}
```

It is not the unrestricted full-space AA_R block.

## Time Evolution

At each RT step:

1. Build or reuse `H_local(I,t)`:

```text
H_local(I,t) = H0_local(I)
             - E(t) . [ R_I delta + xi_local(I) ]
```

2. Exchange only neighbor Wannier coefficients `C_J`.

3. Apply directed face-neighbor blocks:

```text
H_flux(IJ,t) = H_flux0(IJ) - E(t) . xi_flux(IJ)
```

4. Propagate `C`.

The default production propagator should be an exponential action on the local
block plus a controlled neighbor update.  Taylor remains a diagnostic fallback
only.

## Density And Polarization

Density is reconstructed from local-owned Wannier functions and, when needed,
their face-neighbor tails:

```text
rho_I(r) = sum_{a,b in local/neighbor cut} D_ab w_a(r) w_b*(r)
```

Each rank accumulates only the grid region it owns or communicates via the
existing sparse density exchange.  Full `rho(:,:,:)` allreduce across all ranks
is not a production path.

Polarization is computed from the same matrices used in the Hamiltonian:

```text
P = -Tr[D (R_I delta + xi_local + xi_flux)]
```

The observable must not use a different position operator than the propagation.

## Implementation Tasks

### Task 1: Add Formal DG-Wannier Types

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`

Add dedicated arrays for local blocks, neighbor flux blocks, propagated Wannier
coefficients, and neighbor exchange maps.  Do not reuse `global_wannier_position`
as production storage.

### Task 2: Export Or Construct Formal Blocks

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`

Construct:

- `H0_local`
- `xi_local`
- `H_flux`
- `xi_flux`

from the same Wannier gauge.  For the current 2x2x2 validation case, the global
Wannier90 basis may be used to construct the blocks, but the RT file format must
store only local and face-neighbor cut blocks.

### Task 3: Add A Formal Propagator

**Files:**
- Create or modify: `src/rt/dg/rt_dg_integrator_wannier_exp.f90`
- Modify: `src/rt/main_tddft.f90`
- Modify: `src/io/inputoutput.f90`

Add `time_integrator_dg_fragment='wannier_exp'`.

The first implementation may use:

- exact diagonalization/exponential for each local `H_local(I,t)`;
- explicit neighbor flux application using exchanged `C_J`;
- no all-to-all communication.

### Task 4: Add Formal Density/Pz Path

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Modify: `src/rt/dg/rt_dg_fragment_ops.f90`

Use only the formal local/neighbor Wannier blocks.  Remove or demote the
center-only and full-global AA_R diagnostic routes from production selection.

### Task 5: Static Validation Before RT

**Files:**
- Modify or create diagnostic under: `tools/`

Check:

- no-flux LCFO gap;
- projected Wannier `H0_local + H_flux` gap;
- Hermiticity of all blocks;
- norm/unitarity of one zero-field propagation step.

RT is not considered valid if the static gap is already collapsed.

### Task 6: Runtime Verification

Run all commands with `OMP_NUM_THREADS=2`.

Required tests:

```bash
cmake --build build-mpi-eigenexa-wannier-lib -j 4
OMP_NUM_THREADS=2 mpirun -np 8 ./build-mpi-eigenexa-wannier-lib/salmon < short_rt_input
OMP_NUM_THREADS=2 mpirun -np 8 ./build-mpi-eigenexa-wannier-lib/salmon < long_rt_input
```

Acceptance:

- no full global position matrix in propagation;
- no full density allreduce path;
- no all-to-all coefficient exchange;
- `Pz` response does not show the previous artificial very-high-frequency
  zero-cross count;
- `rt iterations / step < 1 s` for the 2x2x2 test.

