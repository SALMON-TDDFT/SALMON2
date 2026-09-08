# Streamed Affine Proof and Point-Group Projection Design

## Goal

Preserve the complete full-system affine symmetry of the global LCFO/MLWF
space without allocating or iterating over a replicated
`Nsym * Nwann * Nwann` representation, while retaining exact operator
projection under the physically relevant point group about a full-system
symmetry center.

## Decision

Use two distinct symmetry layers.

1. The complete affine group, including primitive translations inside the
   supercell, is a proof layer.  SALMON streams point-permuted orbital rows,
   validates identity, metric unitarity, and generator/closure relations, and
   reduces scalar receipts and a canonical fingerprint.  No complete affine
   representation tensor is retained.
2. The maximal subgroup whose operations share a full-system fixed center is
   an operator-projection layer.  Its order is bounded by the crystallographic
   point group (at most 48 in three dimensions).  Only this subgroup may be
   materialized as dense `Nwann * Nwann` matrices for scalar/vector operator
   projection.

Translations are therefore not discarded.  They are proved by the first
layer and omitted only from the dense operator-projection storage where a
common fixed center is required.

## Data flow

The conventional DC calculation and LCFO+EigenExa produce the Gamma-real
global candidate window.  The full atomic structure produces the affine
catalog, point maps, multiplication table, translation subgroup, point
cogroup, and cocycle.  The occupied-plus-buffer seed space is orthonormalized
without materializing all affine images.  Its already-measured full-affine
residual must pass before Wannier90 is entered.

Wannier90 receives only the bounded coordinator M/A matrices and returns the
MLWF transform.  SALMON canonicalizes that transform and then streams the
complete affine action over row-owned MLWF data.  The stream produces identity,
unitarity, closure, and fingerprint receipts.  These receipts bind the full
affine proof to V3.

For observable publication, SALMON derives the maximal exact fixed-center
subgroup directly from the full-system atomic catalog, constructs only its
dense representation, and projects the overlap, Hamiltonian, position, and
velocity matrices.  Position is shifted to the same full-system center during
vector covariance checks.  Fragment boundaries and buffer redistribution do
not define this center.

## Memory and complexity contract

The production path forbids:

- replicated `Nsym * Nwann**2` arrays;
- repeated `Nsym**2 * Ngrid` verification after the affine catalog has already
  established the multiplication table;
- full-wavefunction gathers on the Wannier90 coordinator; and
- fragment-local symmetry centers for full-system operator projection.

Full-affine work is bounded by row/orbital tiles plus point-exchange buffers.
Dense symmetry storage is bounded by `Npoint * Nwann**2`, where
`Npoint <= 48`.  V3 records the full-affine streamed receipts separately from
the fixed-center point-group projection receipts.

## Failure behavior

Publication fails collectively when any of the following occurs:

- the LCFO seed space is not closed under the full affine action;
- the streamed post-MLWF identity, unitarity, or closure receipt exceeds the
  configured symmetry tolerance;
- the affine fingerprint or fixed-center subgroup fingerprint is absent or
  rank-inconsistent;
- no exact inversion exists for a centrosymmetric input such as ideal Si64;
- the fixed-center subgroup representation is nonunitary or fails scalar/vector
  covariance; or
- measured workspace exceeds its declared limit.

There is no fallback to the custom localizer or to a noncanonical checkpoint.
Normal DC LCFO+EigenExa and the generalized-eigenvalue Exp-only V3 RT boundary
remain unchanged.

## Verification

TDD fixtures compare the two-layer path with dense references for small groups
on MPI 1/2/4/8.  Source contracts reject the forbidden all-affine dense tensor
and quadratic action-table/grid proof.  A clean MPI+ScaLAPACK+EigenExa+spglib+
Wannier90 overlay must pass before genuine ideal Si64 is rerun.

Si64 acceptance requires 384 retained and 128 occupied states, full-system
inversion, bounded memory, an accepted V3 checkpoint, field-off stationarity,
impulse linear response, and long-pulse HHG spectra computed from polarization.
Current remains a secondary `dP/dt` consistency check.  H2 and H4 are reported
as peak, dip, or slope and must be suppressed relative to H3 for the ideal
centrosymmetric structure.
