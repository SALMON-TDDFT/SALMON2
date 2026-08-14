# Factored Point-Cogroup Overlap Cache Design

## Problem

The post-Wannier factored point-cogroup proof repeatedly reconstructs the same
row-owned point-representation matrices. Each reconstruction divides 384
orbitals into 32-column tiles and performs one `MPI_Reduce` per destination
rank and tile. With 48 point operations and roughly 188 generator relations,
this creates tens of thousands of small collectives.

## Design

1. Assemble every point-cogroup representation once and retain the row-owned
   tensor for the duration of the proof. Reuse cached left and right matrices
   in every generator relation. Direct expected affine actions remain assembled
   from their spatial maps, preserving the independent cocycle check.
2. Replace the per-owner reductions inside overlap assembly with one
   `MPI_Reduce_scatterv` per orbital tile. Compute the complete local row block
   with a matrix multiplication, then scatter the globally reduced row blocks
   to their owners.
3. Keep the public numerical outputs and proof tolerance unchanged. Account for
   the cached tensor and tile buffers in the checked workspace receipt and use
   collective allocation failure handling.

For Si64 this adds about 14 MiB per rank for 48 cached 384-by-48 row blocks,
while eliminating repeated left/right assembly and reducing inner collectives
by the MPI rank count.

## Verification

- Compare the optimized row-owned overlap with the existing analytical fixture
  on MPI 1, 2, 4, and 8.
- Assert the factored proof prepares exactly 48 point representations rather
  than rebuilding left/right representations for every checked pair.
- Keep the nontrivial cocycle GREEN and corruption RED fixtures.
- Build the production SALMON target and run the route checker.
