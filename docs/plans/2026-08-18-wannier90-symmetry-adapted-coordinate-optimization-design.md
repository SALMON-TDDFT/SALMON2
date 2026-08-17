# Wannier90 Symmetry-Adapted Coordinate Optimization Design

> **Status: superseded on 2026-08-18.** The DMN already contains the complete
> finite-group operation set, so upstream Wannier90's one-pass Reynolds
> average is the exact gradient projector.  The observed cost came from a
> SALMON patch that mistakenly replaced it with repeated averaging.  Commit
> `c08e91f8` restores the simpler upstream algorithm; this decomposition design
> must not be implemented unless a future case proves the complete-group
> projector mathematically insufficient.

## Goal

Minimize the Wannier spread while preserving the supplied finite-group DMN
symmetry throughout optimization.  Replace repeated full-matrix Reynolds
projection with optimization coordinates that contain only symmetry-allowed
rotations.

## Problem

The current patched Wannier90 path computes an unconstrained dense gradient and
then repeatedly averages it with the supplied symmetry generators.  On Si64,
the projection needs more than one hundred dense iterations per Wannier cycle.
The projected gradient is nearly zero and the spread remains near 3413 Angstrom
squared, compared with about 1945 Angstrom squared in the unconstrained Gamma
path.  Increasing iteration caps treats the symptom, scales poorly, and does
not reveal whether the symmetry representation leaves any localization
freedom.

## Mathematical Model

For the finite symmetry group, decompose the supplied Wannier representation
as

```text
D(g) = direct_sum_alpha I(m_alpha) tensor rho_alpha(g).
```

Here `rho_alpha` is an irreducible representation and `m_alpha` is its
multiplicity.  Unitary rotations that preserve the representation are exactly

```text
V = direct_sum_alpha V_alpha tensor I(dim(rho_alpha)),
```

where each `V_alpha` belongs to `U(m_alpha)`.  The anti-Hermitian Lie algebra of
these multiplicity rotations is the complete set of symmetry-allowed Wannier
search directions.

Wannier90 will therefore project a dense spread gradient into multiplicity
coordinates once per cycle, perform conjugate-gradient and line-search work on
those smaller blocks, and expand the accepted update back into the retained
space.  Symmetry is preserved by construction rather than by convergence of an
iterated projector.

## Architecture

At setup time:

1. Read and validate the DMN generator representations.
2. Reconstruct the finite group representation by deterministic generator-word
   traversal without retaining unnecessary duplicate matrices.
3. Separate isotypic components with character projectors.
4. Resolve each component into multiplicity and irreducible coordinates.
5. Validate reconstruction, group products, unitarity, and the band/Wannier
   intertwining relation.
6. Publish multiplicities, allowed Lie-algebra dimension, residuals, memory,
   and a deterministic fingerprint.

At each Wannier cycle:

1. Compute the ordinary spread gradient.
2. Contract it into the allowed multiplicity blocks.
3. Run search-direction and line-search operations on those blocks.
4. Expand the accepted anti-Hermitian update once.
5. Apply the update and measure, but do not iteratively repair, covariance.

No array containing every full group matrix is persistent.  Generator matrices,
isotypic bases, small multiplicity blocks, and one streamed group workspace are
allowed.

## Failure Contract

Fail before localization if:

- the supplied operation catalog is not a finite associative group;
- generator words do not reproduce the catalog;
- any representation is nonunitary or violates a group product;
- the decomposition cannot reconstruct every generator within tolerance;
- band and Wannier representation multiplicities do not admit the required
  intertwiner;
- checked extents, allocation consensus, or MPI metadata agreement fail.

Fail during localization if an expanded update violates covariance or
unitarity beyond tolerance.  Do not increase projection iterations or silently
weaken tolerance.

## Performance Contract

The representation analysis is one-time preprocessing.  Each localization
cycle performs a bounded number of retained-space contractions plus small
multiplicity-block operations.  It must not perform a convergence loop of dense
`N x N` matrix products.  All persistent and transient allocations use checked
wide arithmetic and collective failure handling.

## Staged Delivery

### Phase A: diagnostic only

Compute the decomposition, multiplicities, allowed rotation dimension, and the
norm of the current gradient inside and outside the allowed space.  Do not alter
the production update.

### Phase B: comparison

For one cycle, compare the analytic allowed-space gradient with the legacy
iterated projector.  Require agreement on small exact fixtures and report the
Si64 difference without changing the result.

### Phase C: optimized update

Use multiplicity coordinates for search direction, line search, and update.
Verify spread decrease and covariance together.

### Phase D: simplify

Remove the legacy 100/1000-iteration symmetry projector and its duplicated
state after the new path is proven.

## Verification

- Exact Z2, Z3, Z2 x Z2, and noncommuting S3 representations.
- Multiplicity-one representations with no nontrivial internal rotation.
- Multiplicity-two representations with a known allowed `U(2)` rotation.
- Forbidden gradients project to zero; allowed gradients remain unchanged.
- Rank-disagreeing metadata, malformed products, nonunitary generators,
  allocation failure, and checked-extent adverse cases reject collectively.
- Si64 replay reports the allowed dimension, lowers spread when freedom exists,
  and keeps covariance within tolerance.
- MPI 1/2/4/8 runs return identical fingerprints and numerical receipts.
- Memory receipts prove that no persistent full-group dense tensor is formed.
