# Point-Orbit Center Blocks Design

## Evidence

The corrected complex Jacobi sweep reaches a stationary objective for Si64,
but the resulting frame has point-action leakage `0.54203` at point operation
30, column 11.  The objective is `53.719`, so the projected periodic-position
components are materially noncommuting.  Requiring 288 point-closed Hermitian
matrices to share scalar eigenvectors is therefore stronger than the physical
requirement and cannot be repaired by more sweeps or a looser tolerance.

## Construction

Use the stationary Jacobi frame only to provide deterministic seed vectors.
For each seed column, form its complete point orbit with the compact retained
representation.  Compute the three periodic-position expectation phases of
every orbit vector and cluster those images by periodic center.  Within each
center cluster, construct the numerical span of the orbit images.  The
resulting spans are accepted only when they are mutually orthogonal, their
ranks sum to the full character multiplicity, and every point operation maps
each span wholly into one other span.

Select the first seed, in the existing deterministic column order, that gives
a complete system of center blocks.  Diagonalize the LCFO discriminator only
inside each block.  Exact residual multiplets remain block-internal and are
bound through the existing full-projector fingerprint.

This constructs the required block-monomial point action directly.  It does
not assert that noncommuting projected position components have simultaneous
scalar eigenvectors.

## Failure contract

Reject collectively when no seed produces a complete orthogonal block system,
an orbit-center cluster loses numerical rank, or a point operation leaks
between the constructed blocks.  Report the worst seed, operation, block, and
residual.  Do not relax `dg_ow_symmetry_tolerance` or the final affine proof.

## Memory and MPI

All new arrays are compact retained-space metadata: at most `O(|P| m^2)`
complex values and `O(|P|m)` center/rank metadata.  Inputs are already
replicated and rank-agreed.  Allocation extent arithmetic is checked before
allocation, allocation failure is reduced collectively, and every numerical
branch is reduced before ranks return.

## Acceptance

- A noncommuting periodic-position fixture that defeats scalar joint
  diagonalization produces a complete point-permuted block system.
- Point-frame rotation leaves block projectors and fingerprint invariant.
- Rank loss, incomplete orbit coverage, and corrupt point action reject.
- W90 MPI 1/2/4/8, route checks, Release build, and Si64 MPI 8 pass without
  relaxing the final center-orbit or affine-cocycle gates.
