# Wannier90 Local Gradient Synchronization Design

## Problem

Wannier90 computes the localization search direction in the distributed
`cdq_loc` array.  The symmetry hook currently projects the separate replicated
`cdq` array, but the line search immediately saves and uses `cdq_loc`.  The
projected direction therefore does not participate in the update.

This explains the Si64 checkpoint: the identity covariance defect is at
roundoff while all five nontrivial affine generators have order-one defects.

## Design

Immediately after `internal_search_direction`, gather `cdq_loc` into `cdq`,
broadcast the complete replicated direction, apply
`sitesym_symmetrize_gradient(2, cdq)`, and copy each rank's owned slice back to
`cdq_loc`.  The existing line search then consumes the projected local search
direction without any change to its numerical algorithm.

The synchronization is placed in the Wannier90 external-project patch rather
than SALMON's runtime wrapper because the defect is internal to Wannier90's
distributed optimizer.  This also preserves the correct behavior for more
than one k point; directly projecting only `cdq_loc` would not provide the
complete k-space array required by the symmetry routine.

## Verification

1. Add a source-level patch regression that fails unless gather, broadcast,
   projection, and local-slice restoration appear in that order.
2. Rebuild the patched Wannier90 library and run focused W90 MPI tests on
   1/2/4/8 ranks.
3. Run construction and route regressions.
4. Rerun Si64 through Wannier90 and require the post-Wannier generator
   covariance gate to pass before proceeding.

