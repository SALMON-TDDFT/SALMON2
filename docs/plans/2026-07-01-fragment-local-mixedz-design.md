# Fragment-Local Mixed-Z Propagation Design

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## Goal

Move the DG mixed-Z real-time propagation from the current global mixed split backend to a fragment-local backend while keeping the validated global backend as the reference.

## Current State

The present production path uses `yn_dg_fragment_rt = 'y'` and `yn_dg_mixed_z = 'y'`, but the propagation backend is still `global_mixed_split_backend`. This path has now produced physically plausible Si impulse and laser spectra after the polarization observable was corrected to use the full mixed-Z position operator. It is therefore the best reference implementation for local propagation.

The existing `fragment_local_mixed_split_backend` route already has layout diagnostics and a storage/writeback skeleton, but it is not yet a trusted production propagation path. The field-on local kernel must be developed against the global backend before it is allowed to replace it.

## Chosen Approach: W + Pself + Pneighbor

Use a local mixed block per fragment containing:

- owned Wannier functions `W`
- plane-wave complement on the owner fragment, `Pself`
- plane-wave complement from face-neighbor fragments, `Pneighbor`

This is the accuracy-first route. Earlier Si results showed that dropping mixed `WP/PP` coherence can remove the physical response, so a W-only route is not acceptable as the final target. `W + Pself` is a useful intermediate diagnostic, but the production target should include neighbor P blocks.

## Propagation Data Flow

At each RT step:

1. Keep `global_mixed_split_backend` available as the reference.
2. Build fragment-local blocks with global IDs for W and P slots.
3. Gather local coefficients for occupied propagated states into each fragment-local block.
4. Apply the split propagation locally:
   - field kick: exponentiate `-E(t) * Z_local`
   - field-free phase: apply local flux/mixed Hamiltonian phase
5. Accumulate/write back owner-local W and P coefficients.
6. Compare local coefficients and observables against the global backend in diagnostic mode.
7. Switch production only when field-off and field-on diagnostics are within tolerance.

## Error Handling And Diagnostics

The local backend must not silently fall back. It should return:

- `candidate_available`
- `replacement_applied`
- `bad`
- `replacement_block_reason`

Diagnostics should report:

- local block dimensions: W, Pself, Pneighbor
- field-on coefficient difference vs global backend
- field-off coefficient difference vs global backend
- Pz difference vs global backend
- norm conservation
- optional per-block worst offender

## Acceptance Criteria

The fragment-local backend is considered usable when:

- field-off propagation matches the global backend to numerical tolerance
- field-on short impulse matches the global backend for Pz and coefficient norms
- Si impulse spectrum from local backend agrees with global backend peak positions and order of magnitude
- Si laser HHG using `omega^2 |Pz(omega)|^2` does not introduce new unphysical large even-order peaks beyond the current global reference behavior

## Out Of Scope

- Removing the global backend
- Optimizing communication
- Changing the physics model or Wannier/PW basis construction
- Replacing environment toggles with namelist controls, except where needed for the new backend selector
