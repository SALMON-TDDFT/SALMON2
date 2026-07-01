# Owner-Local Mixed-Z Redefinition Design

## Goal

Redefine the fragment-local mixed-Z propagator so the propagated state is an owner-unique, closed local subspace rather than a block that temporarily includes neighbor P states and discards them at writeback.

## Current Failure Mode

The current fragment-local candidate builds field blocks from `W_owner + P_self + P_neighbor`.  The block is propagated as if all entries are dynamical, but writeback only returns `W_owner` and `P_self` to global storage.  `P_neighbor` amplitudes are used during the field kick and then dropped.

This explains the observed behavior:

- `w_only` has almost no length-gauge response.
- `w_pself` gives a small response with the wrong phase/sign tendency.
- `all` recovers some early-time field correlation but creates a residual polarization shift.
- current global and fragment-local have similar P amplitude scale, but fragment-local loses the clean fundamental response after the pulse.

The block is not a closed time-evolution space, so it should not be treated as a physical propagation backend.

## Chosen Architecture

Use an owner-defined local dynamical subspace:

```text
dynamic state per owner fragment = W_owner + P_owned
```

Here `P_owned` means P channels whose owner fragment is the current fragment.  Neighbor information may enter only as an embedding or matrix-construction ingredient, not as a propagated coefficient that is later discarded.

The writeback rule becomes one-to-one:

```text
W_owner -> global W flux-basis coefficients
P_owned -> global P coefficients
```

No neighbor P coefficient is written by a non-owner fragment.  No propagated coefficient is dropped.

## Owner-Local Result and Halo-Owner Route

The first conservative implementation, `owner_local = W_owner + P_owned`, is a closed writeback space but is too small for the observed linear response.  The short C64 laser check showed that it removes almost all field response, consistent with the diagnostic result that most of the mixed-Z matrix norm is outside the owner-only block.

The production route should follow the original DG Hamiltonian structure:

```text
updated coefficients = owner W + owner P
neighbor coefficients = communicated/read-only halo data
local action = extended Hamiltonian block over W_owner + P_self + P_neighbor
writeback = owner W + owner P only
```

In this definition, neighbor P is not a local owned degree of freedom and is not averaged back.  It is used only as input data when evaluating the owner rows of the local Hamiltonian action.

The explicit input selector for the production route is:

```fortran
dg_mixed_z_frag_local_field_block = 'halo_owner'
```

The older `all` name remains a compatibility alias for the same neighbor-read/owner-write behavior.

## Verification

The required checks are:

- Field-off owner-local run: `Pz_ptp` should remain at numerical noise.
- Laser short run: Pz should correlate with the field/vector-potential proxy at the fundamental.
- Current global vs owner-local short run: compare Pz amplitude and phase for the first 100 au.
- Full laser run: after pulse, residual Pz should be smaller than the in-pulse fundamental response and no artificial DC shift should dominate.
- HHG: H1 must be present before judging H3/H5/H7.

## Non-Goals

- Do not implement ad hoc averaging of neighbor P writeback.
- Do not make `P_neighbor` writeback owner-by-last-writer, sum, or simple average.
- Do not tune the HHG spectrum by changing the observable definition.
