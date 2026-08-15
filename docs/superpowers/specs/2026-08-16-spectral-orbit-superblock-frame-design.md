# Spectral Orbit-Superblock Trial Frame

## Problem

The Si64 spectral basin catalog contains 514 basins in six affine-group
orbits.  The existing rank selector assigns 384 columns in total, but it
diagonalizes each representative basin operator independently and assumes the
propagated columns from different orbits are mutually orthogonal.  The measured
post-propagation Gram defect is exactly `1.0`; equal column count therefore does
not establish a complete orthonormal retained frame.

An arbitrary basin/global-grid numbering must not decide which overlapping
orbit subspace is retained first.

## Physical Orbit Signature

For every affine basin orbit, construct a decomposition-independent signature
from:

1. orbit cardinality;
2. the tolerance-quantized full sorted spectrum of its representative basin
   operator;
3. the tolerance-quantized representative operator trace; and
4. the existing basin/operator provenance fingerprints.

Sort signatures lexicographically.  Signatures equal within tolerance form one
orbit superblock and are processed simultaneously.  Global point IDs and basin
labels may break ties only when emitting receipt names after the numerical
subspace has already been fixed; they must not affect subspace selection.

## Complement Selection

Maintain a distributed accepted frame `Q` whose columns are orthonormal and
whose span is invariant under every supplied affine generator.

For one orbit superblock:

1. stream the full representative eigenspaces needed in descending spectral
   blocks;
2. propagate each complete spectral block through every basin in its affine
   orbit;
3. project the entire superblock candidate matrix against the previously
   accepted frame, `C <- C - Q(Q^H C)`;
4. diagonalize `C^H C`, reject nonfinite values and never split a
   tolerance-degenerate singular-value block;
5. retain complete numerical-rank blocks until the superblock's catalog rank
   is reached; and
6. apply symmetric orthonormalization inside the whole superblock.

Adding a complete affine orbit at once keeps the accepted span invariant.
Symmetric orthonormalization is a function of the covariant Gram matrix, so it
commutes with the induced generator action.  If the requested rank is not
available after projection, reject with the superblock signature, available
rank, and requested rank rather than silently duplicating a direction.

The accumulated accepted rank must finish at exactly the retained-state count
and its global Gram defect must be at most `10*tolerance`.

## Downstream Block Contract

The final columns are grouped by physical orbit superblock, not by individual
basin.  Replace the basin-column offsets used by
`build_dg_spectral_channel_generator_actions` with a `channel_superblock`
catalog (column offsets plus one superblock ID per basin orbit).  Generator and
fixed-center action validation must require zero leakage between distinct
superblocks while allowing dense unitary mixing inside a superblock.

Wannier90 receives the resulting full dense symmetry matrices.  No assumption
that one output column belongs to exactly one spatial basin remains.

## Distribution and Memory

Keep `Q` row-owned as `nlocal x Nstate`.  Process one orbit superblock and one
representative spectral block at a time.  Small Gram, overlap, spectrum, and
polar matrices may be replicated.  Do not retain all basin eigensystems or an
`Npoint x Ngenerator` map.  All extents, LAPACK workspaces, and byte receipts
are checked before allocation and reduced with `MPI_MAX`.

## Collective Failure Contract

Rank agreement is required for signatures, superblock boundaries, requested
ranks, spectra, tolerance, and provenance before shape-dependent collectives.
Allocation, LAPACK, numerical-rank, residual, and MPI failures are reduced
collectively.  Partial outputs are deallocated on every failing path.

## Verification

1. Add a synthetic two-orbit case whose independently selected channels are
   nonorthogonal but whose combined candidate span has full rank.  The old
   propagation must fail its Gram gate; the superblock complement construction
   must return an orthonormal frame.
2. Rotate representative eigenvectors independently inside degenerate spectral
   blocks and require identical final projector/fingerprint.
3. Permute basin numbering and global row ownership and require identical
   superblock signatures, frame projector, and action receipts.
4. Add a tied-signature case and require simultaneous processing; a
   rank-deficient tied superblock must reject collectively.
5. Run MPI 1/2/4/8 construction and EigenExa fixtures, route checks, full build,
   then Si64 through channel propagation and fixed-center action generation.

