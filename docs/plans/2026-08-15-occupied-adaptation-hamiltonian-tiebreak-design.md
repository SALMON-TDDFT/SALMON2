# Hamiltonian-resolved occupied adaptation design

## Problem

The Si8 production run reaches the cocycle-aware point-cogroup occupied
adaptation with stable memory, then rejects because requested rank 16 cuts a
finite degenerate block of the group-averaged occupied projector.  The measured
boundary is

- selected edge `0.31788773991605052`;
- rejected edge `0.31788773991604996`; and
- gap `5.55e-16`.

Changing the tolerance cannot define a physical rank-16 subspace.  Choosing an
arbitrary eigensolver prefix would make the result depend on MPI decomposition,
library details, and the input gauge.

## Decision

Retain the group-averaged occupied projector as the primary discriminator.  If
the requested rank crosses one of its tolerance-degenerate blocks, project a
physical Hamiltonian into exactly that block, average the projected Hamiltonian
over the same validated translation/point-cogroup action, and diagonalize it as
a secondary discriminator.

All primary blocks wholly above the boundary are retained.  The secondary
Hamiltonian supplies only the remaining columns from the boundary block.  If
its own tolerance-degenerate multiplet crosses the requested boundary, reject;
do not use row number, LCFO ordinal, atom centre, spread, or an arbitrary
eigensolver gauge to split it.

The selected subspace remains fixed rank and preserves electron count.  After
selection, remeasure orthonormality, Gamma reality, density drift, and closure
under the complete affine generator set.  Passing the secondary eigensystem is
not a substitute for those final gates.

## Physical Hamiltonian and covariance

The secondary operator must be the Hamiltonian represented in the same smooth,
full-system distributed frame used by the occupied adaptation.  It may be
assembled by streaming Hamiltonian matrix elements or by transporting a
payload-bound upstream Hamiltonian receipt into that frame.  A bare diagonal
array of fragment LCFO eigenvalues is not sufficient after occupied and
point-cogroup rotations.

For an input-frame rotation `Q`, the projector and Hamiltonian transform as
`P -> Q^H P Q` and `H -> Q^H H Q`.  The final selected projector must therefore
be invariant after mapping back to real space.  The Hamiltonian payload,
catalog, cocycle, frame, and selected projector fingerprints are chained into
the resulting receipt.

## Distributed algorithm

1. Build and diagonalize the existing translation-then-point-cogroup averaged
   occupied projector.
2. Identify the complete primary blocks and the unique boundary block using
   collectively agreed adjacent-gap decisions.
3. Stream the physical Hamiltonian into the averaged-orbit coordinates without
   retaining an operation-by-state-by-grid tensor.
4. Project the averaged Hamiltonian into only the boundary block.
5. Diagonalize the small Hermitian block and select complete secondary blocks
   until the requested rank is reached.
6. Reconstruct the selected real-space frame one vector at a time and run all
   existing global receipts and affine closure gates.

The persistent memory remains the existing distributed orbit eigensystem.  New
memory is bounded by one streamed Hamiltonian tile and the boundary-block dense
matrices.  All extents, MPI counts, and byte receipts are checked before
allocation, and allocation or numerical failure is rejected collectively.

## Failure behavior

Collectively reject:

- a non-Hermitian, nonfinite, stale, or rank-disagreeing Hamiltonian;
- a Hamiltonian receipt not bound to the smooth retained frame;
- inconsistent primary or secondary cluster boundaries across ranks;
- a secondary multiplet that is still cut by the requested rank;
- loss of rank, orthonormality, Gamma reality, density conservation, or affine
  closure after reconstruction;
- unsafe extent, MPI count, allocation, or receipt arithmetic.

Diagnostics publish primary selected/rejected edges, primary boundary-block
dimension, secondary selected/rejected edges, residuals, and conservative peak
workspace.  Failure diagnostics are emitted before the collective stop.

## Verification

- A synthetic projector-degenerate block split by a covariant Hamiltonian is a
  GREEN case and returns the same projector after arbitrary input rotation.
- A projector and Hamiltonian that remain jointly degenerate across the target
  boundary are rejected.
- Non-Hermitian, nonfinite, provenance-mismatched, ownership-corrupt, and
  rank-disagreeing inputs are rejected collectively.
- Projector, fingerprint, cluster receipts, and workspace receipts agree on MPI
  1, 2, 4, and 8 ranks.
- Si8/MPI8 passes the currently measured `0.3178877` primary boundary and is
  monitored through the next production stage without relaxing tolerances.

