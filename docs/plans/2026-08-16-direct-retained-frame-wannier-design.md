# Direct Retained-Frame Wannier Design

## Goal

Remove spectral-basin channel selection from the production Wannier route.
Use the already complete, orthonormal retained frame as the single Wannier90
trial frame, while preserving the occupied/complement construction and the
existing symmetry checks.

## Motivation

The Si64 diagnostic run reached the spectral trial Gram gate with every trial
column normalized, but with unit off-diagonal overlap.  In particular, basin
74 column 88 and basin 50 column 25 were the same retained-space direction.
The failure is therefore not missing orthogonalization or insufficient
retained rank.  Independent basin eigensystems and symmetry propagation can
select the same direction more than once, while the sum of their requested
ranks still equals the complement dimension.

The retained frame already consists of:

1. the symmetry-adapted occupied DC-LCFO block; and
2. the complete s/p projector block orthogonalized against that occupied
   block.

It is a complete orthonormal frame before basin processing.  Basin processing
is not needed to make it orthogonal and duplicates work that Wannier90 is
intended to perform.

## Chosen Architecture

1. Preserve construction of the symmetry-adapted occupied block and its
   occupied-orthogonal s/p complement.
2. Use their concatenated distributed retained frame directly as the trial
   coefficient frame.  In retained coordinates this is the identity; on the
   spatial grid the anchors are the retained basis functions themselves.
3. Build the Wannier target symmetry matrices in the same retained frame as
   the band symmetry matrices.  Thus the initial target representation is
   `D_wann = D_band`, with identical operation ordering, orientation, and
   provenance.
4. Invoke Wannier90 once for the complete retained space.  Do not invoke it per
   occupation block, basin, orbit, or character sector.
5. Keep the existing post-Wannier character-sector alignment, Gamma sewing,
   inverse transform, point-cogroup proof, and checkpoint provenance gates.

The equality `D_wann = D_band` is an initial symmetry convention, not a claim
that the output orbitals equal the input retained frame.  Wannier90 may rotate
within the intertwiners allowed by the representation.  Si64 acceptance must
show that this leaves enough freedom to localize; otherwise the fallback is an
algebraically constructed regular target representation, not a return to
independent basin eigensystems.

## Removed Production Work

The production path no longer performs:

- periodic spectral-basin construction;
- preparation and projection of hundreds of basin operators;
- one dense eigensolve per basin orbit;
- independent basin rank selection;
- basin-orbit channel propagation;
- basin-induced target-operation construction; or
- storage of basin labels, spectra, selected ranks, and propagation metadata.

The reusable basin primitives may remain temporarily for focused tests and
historical comparison, but the route checker must prove that production no
longer calls them.

## Data Flow

```text
adapted occupied DC-LCFO states
              +
occupied-orthogonal complete s/p complement
              |
              v
  complete orthonormal retained frame
       |                     |
       |                     +--> D_band == D_wann
       v
spatial Wannier anchors
       |
       v
single symmetry-constrained Wannier90 call
       |
       v
existing character/Gamma/inverse/cocycle validation
```

No replicated real-space `N x M` trial frame or all-basin operator collection
is introduced.  Trial coefficients are generated analytically as row-owned
identity rows or bypassed where the existing materialization API can consume
the retained frame directly.

## Invariants and Failure Handling

- The retained frame dimension is `ntarget` on every rank.
- Global retained-row IDs are owned exactly once.
- The retained Gram defect and occupied/complement cross-Gram defect pass
  before Wannier90 preparation.
- `D_band` and `D_wann` use the same canonical operation ordering and action
  direction.
- The target representation fingerprint is bound to the retained-frame and
  operation-catalog fingerprints.
- Every allocation and MPI count added by this route uses checked wide
  arithmetic and collective allocation consensus.
- Trivial translation groups and `ntarget == nstate` remain valid.
- Failure before Wannier90 leaves no partially initialized Wannier workspace.

## Performance and Memory

For real-space size `N`, retained dimension `M`, basin orbit count `B`, and
MPI rank count `P`, this removes approximately

`O(B*N*(M-nstate)^2/P + B*(M-nstate)^3)`

preprocessing work, as well as basin catalog and eigensystem storage.  The
remaining leading work is the single Wannier90 solve and existing distributed
symmetry processing.  Per-rank trial-frame memory remains `O(N*M/P)` through
the retained frame already required by production; no second full frame is
kept.

## Test Strategy

1. Add a RED route test requiring production to avoid every spectral-basin
   constructor, projector, selector, and propagator.
2. Add a focused direct-frame fixture that verifies:
   - row-owned identity coefficients;
   - retained Gram preservation;
   - exact `D_wann = D_band` orientation;
   - canonical fingerprints across MPI 1/2/4/8; and
   - a nontrivial symmetry-allowed unitary rotation within a repeated irrep.
3. Preserve the diagnostic duplicate-basin fixture as evidence that the old
   algorithm fails, but do not run basin processing in production.
4. Run construction, EigenExa, W90, route, obsolete-route, and overlay build
   checks.
5. Run Si64 through Wannier90 completion.  Acceptance requires finite spread,
   convergence or an explicit iteration-limit diagnostic, retained unitarity,
   post-gauge symmetry receipts, and stable per-rank memory.
6. If Si64 is symmetry-overconstrained with `D_wann = D_band`, stop and design
   an algebraic regular target representation.  Do not silently weaken or
   remove symmetry constraints.

