# DG+PW Orthogonalized Mixed-Basis Design

## Goal

Make `DG+PW` follow a single, consistent representation:

- always use the orthogonalized mixed basis produced by `diagonalize_mixed_basis`
- remove the raw fragment-plane-wave propagation path from the runtime hot path
- evaluate fragment-plane-wave off-diagonal couplings from `core+buffer` support integrals

This keeps the runtime representation aligned with the documented design intent in `ADAPTIVE_BASIS_IMPLEMENTATION.md` and avoids mixing raw nonorthogonal coupling with orthogonalized propagation.

## Context

The current codebase has two distinct `DG+PW` paths:

1. A raw coupling path that builds `S_mat_frag_pw` and `H_mat_frag_pw`, keeps fragment coefficients as-is, starts `coef_pw=0`, and lets PW components grow during propagation.
2. A diagonalized mixed-basis path that forms the generalized eigenproblem, applies Lowdin orthogonalization, and stores the resulting transform in `mixed_transform` / `coef_mix`.

These are both present in runtime code today. That means current calculation, density reconstruction, energy evaluation, and propagation can operate on partially different mathematical objects depending on stage and options.

## Design Decision

Adopt **orthogonalized mixed-basis only** for `DG+PW`.

That means:

- `DG+PW` startup must always build and diagonalize the mixed basis before propagation.
- Runtime propagation, current, density, and energy paths must use the mixed-basis coefficients as the authoritative state.
- Raw fragment/PW coefficient propagation remains available only as an implementation detail during startup or internal reconstruction, not as an alternate runtime physics path.

## FP Coupling Design

Fragment-plane-wave off-diagonal matrix elements should be evaluated from the fragment real-space support used for the basis itself:

- integrate on `core+buffer` support for fragmented axes
- use fragment-box periodic wrapping on non-fragmented axes
- do not expand support with fake buffer on axes where `num_fragment(dir)=1`

This matches the axis-wise buffer interpretation already established for buffered DG work:

- fragmented directions carry physical fragment buffer
- non-fragmented directions are plain periodic directions of the fragment box

## Data Model

For `DG+PW`, the authoritative runtime state becomes:

- `mixed_transform`
- `coef_mix`
- `mixed_basis_dim`
- `mixed_basis_ready`

Raw storage remains only for support:

- fragment basis functions in real space
- PW basis metadata
- raw overlap / Hamiltonian couplings used to build the mixed problem

But once `mixed_basis_ready=.true.`, runtime operators should derive their action from the mixed basis, not from the raw fragment/PW split.

## Runtime Rules

### Startup

- If `use_plane_wave_basis=.true.` and `n_plane_waves>0`, call `diagonalize_mixed_basis` at startup.
- Do not use `prepare_mixed_basis_startup` as a steady-state fallback for propagation.

### Propagation

- Propagate the mixed coefficients.
- Any raw fragment/PW arrays needed by lower-level code must be synchronized from the mixed state in a single, well-defined place.

### Density and Observables

- Current calculation should prefer mixed-basis operator application.
- Density reconstruction should use the mixed basis when `mixed_basis_ready`.
- Energy evaluation should not fall back to raw FP coupling when mixed-basis mode is active.

### FP Off-Diagonal Integrals

- Recompute fragment-PW overlap and Hamiltonian coupling using `core+buffer` support.
- Ensure axis-wise support rules are respected:
  - fragmented axis: `core+buffer`
  - non-fragmented axis: local periodic box

## Error Handling

Fail fast for inconsistent `DG+PW` setup:

- `use_plane_wave_basis=.true.` but `mixed_basis_ready=.false.` after startup
- runtime path entering raw fragment/PW propagation when mixed-basis mode is required
- shape mismatches between mixed-basis dimension and fragment/PW backing arrays

Warnings are not sufficient here because silent fallback recreates the representation-mixing bug we want to remove.

## Testing Strategy

### Functional

- `DG+PW` startup reaches mixed-basis diagonalization and sets `mixed_basis_ready=.true.`
- runtime does not use raw fragment/PW current path when mixed-basis mode is active
- density reconstruction works with mixed basis only

### Physics Sanity

- `DG+PW` current for `x/y/z` polarization is finite and stable
- current anisotropy is not dominated by representation mismatch
- `A(t)` is identical between Full and DG runs at matched times

### Regression

- pure `DG` remains unchanged
- `DG+PW` with zero plane waves behaves like pure `DG`

## Implementation Notes

- This is a representation cleanup, not a request to solve all remaining superlattice physics issues.
- For solids, this design is closer to the original intended use: orthogonalized hybrid basis for metallic or delocalized response.
- The immediate success criterion is consistency of runtime representation, not perfect agreement on every one-body energy term in superlattice stress tests.

## DG Hamiltonian Eq.(4) Mapping Note

The SC'25 DG-TDDFT paper describes the DG Hamiltonian in Eq.(4) as a sum of:

1. volume kinetic term
2. volume local-potential term
3. nonlocal pseudopotential term
4. DG face terms (average/jump and penalty terms)

The current SALMON implementation matches that structure only partially.

### What is already consistent

- The block sparsity pattern uses `3^Dim` nearest-neighbor adjacency.
  In 3D, this means up to 27 row-block couplings including self.
- Diagonal fragment blocks are built from volume kinetic and local-potential integrals.
- Off-diagonal DG surface coupling is treated only for face-neighbor pairs, which is consistent with the paper's statement that the boundary integral terms affect edge-adjacent blocks only in the block graph sense of DG interfaces.

### What is not identical to the paper-level Hamiltonian assembly

- The nonlocal PP contribution is not fully assembled into `H_mat` as a persistent DG Hamiltonian block matrix.
  Instead, it is applied later during propagation / observable evaluation through projector application paths.
- Therefore, the runtime one-body operator is mathematically closer to
  `H_local_blocks + H_face_blocks + on-the-fly nonlocal application`
  than to a single monolithic Eq.(4) matrix assembled once.

### Important code-status note

`src/rt/dg/rt_dg_fragment_hamiltonian.f90` still contains legacy comments around the missing/next-stage face-term assembly.
The present code does inject SIPG-like face terms for face-neighbor halo pairs, but the comments are more conservative than the current implementation state.
This should be cleaned up later so that code comments and actual assembly path say the same thing.

### Practical interpretation

- For discussions about 27-neighbor communication/block structure, the current implementation is broadly aligned with the paper.
- For discussions about whether the Hamiltonian matrix is assembled exactly as Eq.(4), the answer is still "not completely":
  the nonlocal term is handled outside the stored DG Hamiltonian blocks.

This distinction matters when comparing energy definitions, current evaluation, and block-level diagnostics against the paper description.
