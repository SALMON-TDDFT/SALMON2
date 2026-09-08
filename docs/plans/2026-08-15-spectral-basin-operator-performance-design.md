# Spectral basin operator performance design

## Goal

Reduce the arithmetic overhead and MPI volume of the streamed spectral-basin
operator without increasing asymptotic memory or changing its numerical
contract.

## Decision

Replace the point-by-point complex outer-product loops with fixed-size packed
tiles and Hermitian rank-k updates.  Each tile stores
`sqrt(point_weight) * state_values(:,point)` and contributes through `ZHERK`.
The tile has a compile-time bounded width, so workspace remains `O(Nstate)` in
addition to the reusable dense basin operator.

Only the upper triangle is accumulated and packed.  MPI reduces the
`Nstate*(Nstate+1)/2` complex values, after which every rank reconstructs the
lower triangle by conjugation.  This nearly halves collective payload relative
to the current dense `Nstate**2` reduction.

Alternatives rejected:

- one `ZHER` call per point retains excessive BLAS-call overhead;
- packing every point in a basin gives good BLAS efficiency but makes peak
  memory depend on the largest basin;
- reducing a full dense matrix preserves unnecessary lower-triangle traffic.

## Contracts

- Prepared point indices remain the only basin traversal metadata.
- Zero weights produce zero packed columns without division.
- Tile and packed-triangle extents are checked before allocation.
- Allocation failure is reduced collectively and cleaned on every rank.
- The optimized result must agree with the standalone reference within the
  existing tolerance on MPI 1/2/4/8.
- The workspace receipt includes the tile and packed triangle.
- A communication receipt reports exactly `Nstate*(Nstate+1)/2` complex
  elements per basin.

