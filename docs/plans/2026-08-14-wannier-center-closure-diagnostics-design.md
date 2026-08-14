# Wannier Center Closure Diagnostics Design

## Goal

Distinguish a genuine post-W90 affine-gauge failure from an ill-defined periodic
center caused by diffuse Wannier functions, without weakening acceptance.

## Design

Keep the exact bipartite center-orbit matching and its existing tolerance.  When
matching fails, compute the failed mapped center's nearest periodic Chebyshev
distance to any returned center.  Include the affine operation index, source
center index, nearest residual, configured tolerance, and the minimum periodic
moment magnitude of that source orbital in the error message.

Pass the already-computed periodic moment magnitudes from production into the
validator as an optional diagnostic payload.  Publish the global minimum and
maximum magnitude before validation so a failed run still records whether its
centers are well-defined.

## Verification

Extend the existing broken-orbit RED to require the new diagnostic fields, then
run construction MPI 1/2/4/8, the route checker, Release build, and diff check.
Acceptance behavior and center ownership remain unchanged.

