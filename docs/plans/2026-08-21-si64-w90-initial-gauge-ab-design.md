# Si64 Wannier90 Initial-Gauge A/B Design

## Goal

Determine whether the Si64 symmetry-adapted Wannier90 stagnation is caused by
the deterministic spectral initial gauge by comparing it with the previously
successful random-projection initialization.

## Fixed conditions

Both runs retain the same SALMON DC/LCFO state, 384-band retained space,
Gamma-only `.dmn` symmetry payload, `site_symmetry=.true.`,
`symmetrize_eps=1d-10`, 200-iteration limit, MPI decomposition, and strict
Wannier90 convergence-receipt gate.  Only the initial projection mode changes.

## Interface

Add a narrow DG input setting `dg_ow_w90_initial_projection` with values
`spectral` (default) and `random`.  The setup routine writes no projections
block for `spectral`, and writes the standard Wannier90 `random` projections
block for `random`.  The generated `.win` is the direct experiment receipt.

## Acceptance

Unit and route tests must prove both generated forms and reject an unknown
mode.  The Si64 A/B input selects `random`; the production run is accepted only
if Wannier90 itself emits its convergence-satisfied receipt.  No Hybrid-SCF
conclusion is drawn from an iteration-limit result.

