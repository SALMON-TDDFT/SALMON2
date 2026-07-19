# WPW Density-Carrying Seed and Fixed-H Relaxation Design

## Status and decision

This document supersedes the occupied-subspace initialization part of
`2026-07-19-wpw-fixed-hamiltonian-basis-flux-relaxation-design.md`.

The first WPW checkpoint keeps the converged DC density and every non-interface
Hamiltonian contribution fixed. It does not run a WPW density SCF. The initial
subspace is built from the occupied fragment orbitals that carry the converged
DC density, followed by a true projection into the nonorthogonal W+P basis.

It must not be described as a projection of the global LCFO eigenvectors.
Those vectors (`coef_wf`) are currently produced only by the Flux-inclusive
`diag_eigenexa` route that this work is intended to avoid.

## Why this source is used

The converged DC state is an ensemble of occupied fragment orbitals. Their
piecewise core-domain density is the density that is frozen for the first WPW
checkpoint. Using their direct sum as the seed therefore preserves the intended
density provenance without first solving the expensive Flux Hamiltonian.

The alternatives were rejected:

1. Running `diag_eigenexa` first would reintroduce the slow Flux-inclusive
   global eigenproblem and use a state already affected by the operator being
   introduced.
2. Treating W overlaps directly as coefficients and setting P coefficients to
   zero is not a projection in a nonorthogonal W+P basis.
3. Reconstructing global orbitals from density alone is not unique and would
   introduce an uncontrolled gauge choice.

## Mathematical contract

Let `B=(W,P)` be the distributed W+P basis and `S=B^dagger B`. For every
occupied fragment orbital `phi_a`, construct both overlap blocks

```text
b_W = W^dagger phi_a
b_P = P^dagger phi_a
b   = (b_W,b_P).
```

The projected coefficients are the rank-qualified solution of

```text
S c = b.
```

The implementation must not replace this solve with coefficient
normalization. After the metric solve, the occupied columns are
S-orthonormalized with rank revelation. Loss of occupied rank is fatal.

Projection diagnostics are:

- relative metric-equation residual `||S C-B||/||B||`;
- captured norm `Tr(B^dagger C)` relative to the source core norm;
- occupied effective rank;
- S-orthogonality defect;
- occupation-derived charge;
- S-metric projector change under an occupied unitary rotation.

## Extra-state contract

Only occupied columns come from the DC density-carrying ensemble. Extra states
are deterministic. Before normalization, every extra column is projected out
of the occupied space:

```text
Q_extra <- Q_extra - Q_occ (Q_occ^dagger S Q_extra).
```

The extra block is then rank-qualified and S-normalized. Random replacement of
an occupied column is forbidden.

## Fixed operator and continuation

The fixed generalized eigenproblem is

```text
(H0_DC + lambda_interface H_interface) C = S C E.
```

Frozen H0 contains:

- WW kinetic, potential, and nonlocal blocks;
- WP volume and nonlocal blocks;
- PP volume and nonlocal blocks.

The complete interface block contains only:

- `ww_face_self`;
- `ww_cross_value`;
- `wp_h_face`.

Support IDs, ownership maps, halo records, and routing schedules are transport
metadata and are never scaled.

The first solve uses `lambda_interface=0`. Accepted solutions are continued to
one with rollback and step reduction. Density reconstruction is diagnostic
only and must never call `wpw_potential_step` in fixed-H mode.

## Frozen-state invariant

The snapshot must cover values, not only stored fingerprints:

- total density, Hartree, XC, and local ionic potential arrays;
- all WW component arrays and IDs;
- WP/PP volume, nonlocal, and face arrays;
- bounded H0 and interface caches;
- layout/ownership IDs and operator provenance.

Validation must be collective and shape-safe. No rank may return before a
collective entered by its peers. Allocation failure must invalidate the whole
candidate and leave no publishable checkpoint.

## Failure and publication policy

Stop before checkpoint publication if any of the following occurs:

- source occupation count or charge is inconsistent with DC provenance;
- W or P overlap construction is incomplete or nonfinite;
- the metric solve fails or loses occupied rank;
- projection residual or S-orthogonality exceeds the active tolerance profile;
- any frozen value, cache, ID, or provenance record changes;
- a fixed-H solve or continuation step is nonfinite or unconverged;
- `lambda_interface` has not reached one.

Long-time HHG and spectral interpretation remain out of scope.
