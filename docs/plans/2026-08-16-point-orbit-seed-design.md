# Point-Orbit Seed Preservation Design

## Context

After Wannier90 converged at iteration 160, Si64 reached the joint periodic-center
gauge and failed to construct complete point-orbit blocks.  In
`build_point_orbit_blocks`, the point-operation loop both constructs orbit
vectors and mutates `orbit_residual` during orthogonalization.  Consequently,
operation `p+1` acts on the residual left by operation `p`, rather than every
operation acting on the same canonical seed.

## Decision

Split the point loop into two stages.  First fill every `orbit_vectors(:,p)` by
applying `point_representations(:,:,p)` to the unchanged normalized seed.
Second, process those stored vectors for center clustering, within-cluster
orthogonalization, and global cover construction.  Reuse the existing
`orbit_vectors` allocation and add no persistent or transient workspace.

## Verification

Add a three-operation cyclic point representation fixture that requires every
operation to act directly on the same seed.  Verify the old implementation
rejects or produces the wrong orbit and the two-stage implementation passes.
Run the Wannier90 MPI fixture on 1, 2, 4, and 8 ranks, route checks, production
build, and then rerun Si64 through the former periodic-center failure.

