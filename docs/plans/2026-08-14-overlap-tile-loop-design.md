# Overlap Tile Loop Optimization Design

## Goal

Reduce the post-Wannier factored symmetry-overlap runtime without increasing its
large-system memory order or changing numerical results.

## Findings

The hot overlap loop contains no floating-point division and its explicit loops
already traverse the first Fortran array index contiguously.  The avoidable work
is instead:

- allocating and freeing four tile arrays for every symmetry and orbital tile;
- materializing `local_tile`, then transposing all of it into `packed_tile`;
- rebuilding full-tile receive counts for every tile;
- issuing more collectives than necessary when the tile width is too small.

The generator-mask branch is bounded by the at-most-48-element point group and
is not a material inner-loop cost.

## Design

Keep the row-owned `MPI_Reduce_scatter` algorithm and its numerical ordering.
Allocate maximum-width tile buffers once per overlap assembly call and reuse
their active slices for every symmetry operation and tile.  Compute the BLAS
product directly in the tile-major layout required by `MPI_Reduce_scatter`,
eliminating `local_tile` and its transpose copy.

Use a compile-time tile width of 64 initially.  This halves collective count
relative to width 32 while adding only bounded tile workspace, not an
`O(N^2)` or all-sector allocation.  Preserve a distinct active width for the
final partial tile and precompute the full-width receive counts.

## Verification

- Add a source-structure regression that rejects per-tile allocation and the
  obsolete `local_tile` transpose path.
- Run the construction MPI fixture on 1, 2, 4, and 8 ranks.
- Run the production route checker.
- Rebuild the Release binary and run `git diff --check`.
- Measure the optimized Si64 run separately because an already-running process
  retains the old executable mapping.

