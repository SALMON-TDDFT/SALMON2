# WPW canonical-face trace assembly design

## Scope

This checkpoint connects real-grid Wannier and windowed plane-wave traces to
the existing SIPG point evaluator.  It assembles only the WP face contribution
on canonical faces.  PP remains volume-only because each windowed plane wave
is a periodic H1 function.  Nonlocal pseudopotential assembly is out of scope.

## Architecture

The face scanner owns geometry and accumulation, but not trace storage or
communication.  A bound trace-provider callback supplies values and gradients
for one face point and both sides.  This keeps rank-local cache layout and the
future support-halo exchange outside the mathematical quadrature routine.

The scanner receives a canonical face identity `(K-, K+, axis, side)`, derives
the outward canonical normal and the unwrapped face grid, and visits every
tangential grid point exactly once.  At each point it requests W and P traces,
checks the provider status, and calls `assemble_wpw_canonical_face_point` with
`h_normal=hgs(axis)`, `sigma=10/h_normal`, and the tangential area weight.
The resulting WP matrix is accumulated into a caller-owned face block.

No global fragment loop, global dense matrix, owner table, or hidden MPI call
is permitted in the scanner.  A remote trace unavailable on the current rank
is an explicit provider failure; the later halo/provider layer must make that
trace available before assembly.

## Geometry and uniqueness

Only entries produced by `build_wpw_canonical_face_list` are admissible, so
`K- < K+` and each physical face is assembled once.  The two wrapped faces of
a two-fragment periodic direction remain distinct through
`side_from_k_minus=-1,+1`.  Grid coordinates remain unwrapped until the point
adapter/provider maps them into periodic storage.  The normal is zero except
for `normal(axis)=side_from_k_minus`.

For a face normal to `axis`, the two tangential ranges are the intersection of
the fragment face extents.  A missing or inconsistent intersection fails
closed rather than silently dropping quadrature points.

## Provider contract

The provider returns, for fixed sorted W row ids and sorted `(K,G)` column ids:

- `W-`, `W+`, `grad W-`, and `grad W+`;
- `P-`, `P+`, `grad P-`, and `grad P+`;
- a nonzero status for missing ownership/halo data, stale cache data, invalid
  ids, or non-finite values.

The scanner validates all result extents and finite values.  The existing face
point evaluator additionally enforces the H1 P-trace equality.  The callback
binding uses an explicit context/lifetime object; no procedure-local target may
escape its lifetime.

## Failure behavior

Invalid face identity, axis, side, grid extent, spacing, weight, callback
binding, callback result, or non-finite trace fails before candidate packing.
No partial face block is returned: output is reset to zero on every failure.
MPI-fatal policy belongs to the outer distributed assembly layer, which can
report the detecting rank before collective termination.

## Tests

1. A two-fragment periodic numerical fixture proves both canonical faces are
   visited once and remain distinct.
2. A deterministic provider proves tangential point count, normal, grid
   coordinates, weight, and accumulated SIPG value.
3. Orientation reversal of traces and normal leaves the face integral
   invariant.
4. Missing provider data, non-finite data, invalid geometry, and P-trace
   mismatch fail closed with a zero result.
5. Source-contract checks reject PP face outputs, global fragment scans, dense
   matrices, and hidden MPI calls.

