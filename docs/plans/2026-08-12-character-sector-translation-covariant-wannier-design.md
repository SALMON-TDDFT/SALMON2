# Character-Sector Translation-Covariant Wannier Design

## Status

This design supersedes the center-matching construction in
`2026-08-12-translation-covariant-post-mlwf-gauge-design.md`.  Center matching
cannot define orbital correspondence when several Wannier functions share a
center, as is normal for complete s/p/d shells.  Translation correspondence
must be derived from the wavefunction representation.

## Physical requirement

The ground-state basis is a 384-dimensional global LCFO space closed under the
actual 32-element atomic translation subgroup and the 48 point-cogroup
representatives.  Wannier90 constrained only by the fixed-center subgroup gives
a useful localized reference gauge, but not a translation-covariant final
gauge.  Translation symmetry is imposed on the fixed basis before RT; RT state
coefficients remain free to break it under a spatially varying field.

## Representation construction

The atomic translation subgroup is finite and abelian.  Stream the translation
action in the accepted LCFO space and retain only a minimal generator set.  The
generators are commuting unitary matrices.  Validate identity, inverse,
commutator, order, and product relations before using them.

Construct the complete character table directly from the validated finite
group product table.  Do not assume a particular decomposition such as
`Z2^5`; arbitrary finite abelian translation groups must work.  Characters,
group elements, and generators receive canonical, rank-independent ordering
and fingerprints.

## Character-sector decomposition

Simultaneously reduce the commuting translation generators.  Each of the 32
characters must occur with multiplicity 12 for the Si64 rank-384 acceptance
case.  More generally, every retained character must have the same internal
multiplicity; otherwise a complete translation orbit basis cannot be formed at
the requested rank.

Use character projectors or an equivalent generator eigenspace refinement.
Process one character at a time and distribute its vectors in the existing
orbital/spatial layout.  Never retain `Ntranslation*Nwann^2` matrices or all
character projectors.

## Wannier90-referenced sector gauge

Wannier90's output supplies the localization reference, not the symmetry
definition.  Project the MLWF reference vectors into every character sector.
Choose a deterministic reference sector and align the other 12-dimensional
sector frames to it using cross-sector localization link matrices and unitary
polar factors.  Degenerate singular-value clusters are treated as complete
blocks; rank loss or a split cluster is rejected.

The alignment must be invariant under input orbital numbering and arbitrary
unitary rotations inside the pre-aligned sector frames.  Gamma-real input must
produce conjugate-paired sectors and a real final localized basis to numerical
tolerance.

## Inverse finite Fourier transform

For each of the 12 aligned internal channels, apply the inverse character
transform over the 32 sectors.  This produces 12 representatives and their 32
exact translates.  The translation action is then a known permutation of
orbit members, including phase conventions fixed by the canonical character
table.  Multiple orbitals may share a center without ambiguity because their
12-dimensional internal channel is resolved before localization.

Apply the resulting unitary transform identically to core and buffer values,
gradients, operator transformations, occupations, and checkpoint provenance.

## Point-cogroup covariance

Each point-cogroup representative permutes translation characters.  Validate
this permutation and its internal 12-dimensional action using the existing
translation cocycle.  Check representative products as point operation plus
translation.  This proves the full 1536 affine closure without materializing
1536 dense matrices.

## Memory model

Allowed persistent data are distributed basis blocks, a small translation
catalog, `O(Nwann)` ownership/sector metadata, minimal generator matrices, and
one character-sector workspace.  One dense `Nwann^2` operation may be gathered
on the coordinator at a time.  Disallowed data include all translation
representations, all affine representations, all character projectors, or all
translated real-space orbitals.

All allocation extents, BLAS dimensions, and MPI counts are checked in a wide
integer kind.  Peak persistent and transient bytes are reported separately.

## Failure handling

Reject nonabelian or nonfaithful translation catalogs, inconsistent generator
orders, incomplete character tables, unequal sector multiplicities, nonunitary
actions, singular localization links, split degeneracy clusters, Gamma reality
loss, density drift, rank change, or cocycle failure.  No symmetry or center
tolerance is relaxed to pass the genuine system.

## Acceptance

- Finite-group fixtures include `Z2`, `Z2 x Z2`, `Z4`, and an adverse
  nonabelian/corrupt table.
- Complex gauge rotations and repeated centers do not change the canonical
  localized result or fingerprint.
- Translation identity, unitarity, commutator, character orthogonality, sector
  multiplicity, inverse-transform, and permutation residuals pass on MPI
  1/2/4/8.
- Density, electron count, rank, and buffer-supported gradients are preserved.
- The 48 representative plus cocycle proof recovers full affine closure.
- The ideal undisplaced Si64 run passes fixed-center DMN, Wannier90,
  character-sector relocalization, center orbit, fragment redistribution, and
  V3 checkpoint.
- Normal DC LCFO+EigenExa and Exp-only RT routes remain unchanged.

Only after V3 acceptance may polarization-primary linear response and
long-pulse Exp HHG resume.
