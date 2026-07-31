# DG PW-Augmented Basis Design

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## Goal

Improve the DG-RT paramagnetic response by augmenting the DC-LCFO-Flux basis with fragment-compatible plane-wave helper functions, instead of switching to real-space propagation.

The target physics is the diamond impulse response:

- `J_tot,z` starts from the kick.
- The paramagnetic current cancels the diamagnetic contribution.
- `J_tot,z` then oscillates around zero.
- Persistent oscillation around `J_dia` indicates missing paramagnetic response.

## Decision

Use the coefficient-space DG formulation and improve the basis. Do not pursue real-space RT propagation for this path.

The real-space route creates a boundary consistency problem: finite-difference stencil and nonlocal pseudopotential application need buffer or ghost wavefunction values, but those buffer values are not independent propagated degrees of freedom. In a proper DG method, ghost values are refreshed from the current stage state and inter-element coupling enters through numerical fluxes. Making real-space RT consistent would therefore require a new DG-FD boundary treatment and PP halo design.

The PW-augmented DG route keeps the Hamiltonian, overlap, momentum/current, and time propagation in the same Galerkin coefficient space. This is closer to the existing DG-RT implementation and targets the likely failure mode directly: the DC-local basis is missing high-energy response channels needed for the f-sum and `J_para`.

## Architecture

The propagated basis is:

```text
B = { DC-LCFO-Flux fragment basis } + { fragment-compatible PW helper basis }
```

The RT state remains coefficient based:

```text
psi(t) = sum_i C_i(t) B_i
```

The time evolution remains:

```text
C(t+dt) = exp(-i S^{-1} H dt) C(t)
```

or the existing Taylor/RK variants in the generalized DG basis.

The Hamiltonian must stay fixed within each Taylor operation. PW augmentation changes the basis and matrix representation, not the meaning of `H` during a Taylor step.

## Existing Code To Reuse

The code already contains several PW augmentation hooks:

- `src/io/salmon_global.f90`
  - `yn_plane_wave_basis`
  - `n_plane_waves_dg`
  - `k_cutoff_plane_wave`
  - `yn_dg_subspace_diag`
  - `dg_subspace_extra_states`
  - `dg_subspace_pw_vectors`

- `src/rt/dg/rt_dg_plane_wave.f90`
  - PW selection and initialization
  - fragment/DC-PW overlap blocks
  - PW-PW overlap and Hamiltonian helper blocks
  - basis compaction support

- `src/rt/dg/rt_dg_fragment_ops.f90`
  - coefficient operations already account for `coef_pw`

- `src/rt/dg/rt_dg_integrator_taylor.f90`
  - Taylor propagation already branches when PW coefficients are present.

The implementation should first audit these paths and fix correctness gaps before adding new machinery.

## Required Cleanup

The abandoned DG Flux real-space RT implementation must be removed before continuing. Revert the commits that introduced:

- `yn_dg_flux_realspace_rt`
- `dc_lcfo_flux_operator.f90`
- complex real-space Flux `hpsi` support
- real-space RT design and plan docs

This keeps the branch conceptually clean and avoids two competing RT designs.

## Basis Construction

The PW helper basis should be fragment-compatible:

- PW candidates are generated from reciprocal vectors up to `k_cutoff_plane_wave`.
- A deterministic subset of at most `n_plane_waves_dg` is selected.
- PW helper functions are projected into fragment-local representation.
- DC-LCFO and PW helper blocks build a generalized overlap matrix `S`.
- The PW helper space should be compacted if it is linearly dependent on the DC basis.

The first production target is not a large PW basis. It is a small sweep that reveals whether extra response channels move the impulse current in the right direction.

## Matrix Requirements

For the augmented basis, these matrices must be consistent:

- `S_DC,DC`, `S_DC,PW`, `S_PW,PW`
- `H_DC,DC`, `H_DC,PW`, `H_PW,PW`
- momentum/current blocks for observables
- nonlocal PP contributions for mixed and PW blocks
- Flux surface terms for all blocks that touch inter-fragment coupling

The same overlap convention must be used in propagation, density reconstruction, observables, and unitarity diagnostics.

## Validation

Use a short diamond impulse case first.

1. Baseline with `yn_plane_wave_basis='n'`.
2. PW helper runs with small `n_plane_waves_dg`, for example 8, 16, 32.
3. Track:
   - `n_plane_waves`
   - overlap condition / compacted PW count
   - RT f-sum proxy or transition strength trace if available
   - `J_para,z`, `J_dia,z`, `J_tot,z`
   - `J_tot,x/y` leakage
4. Success criterion:
   - `J_para,z` cancels more of `J_dia,z`
   - `J_tot,z` moves toward zero-centered oscillation
   - no new large transverse current appears
   - norm/unitarity diagnostics do not degrade

## Risks

- Existing PW augmentation may currently be wired for experiments, not production.
- Mixed `H_DC,PW` and `P_DC,PW` blocks may be incomplete or inconsistent with Flux.
- PW basis can become nearly linearly dependent on the DC basis; compaction and overlap conditioning are required.
- Increasing PW count may improve f-sum but destabilize propagation if `S^{-1}H` is ill-conditioned.
- Nonlocal PP mixed blocks are likely the highest-risk matrix component.
