# Hybrid Hpsi Spatial Redistribution Design

## Problem

The row-owned overlapping-Wannier grid IDs (`ow_core_ids`) are not distributed
by the same MPI ownership as the total-system real-space grid (`dc%mg_tot`).  The
current full-cell Hamiltonian callback incorrectly assumes these layouts are
identical and rejects valid Si64 data before calling `hpsi`.

## Decision

Build and cache a bidirectional MPI redistribution schedule once, then reuse it
for every bounded orbital tile and every hybrid-SCF Hamiltonian rebuild.  Do not
materialize a replicated full-cell orbital array.

## Data Flow

For each tile of at most 16 orbitals:

1. Pack values in local `ow_core_ids` order.
2. Use the cached forward `MPI_Alltoallv` schedule to place values in local
   `dc%mg_tot` ownership order.
3. Apply SALMON's established total-system `hpsi` operator.
4. Use the cached reverse schedule to return `H psi` values to the original
   `ow_core_ids` owners and ordering.
5. Accumulate the distributed projected Hamiltonian rows.

The schedule is keyed by the communicator, global grid extent, source IDs, and
destination IDs.  It is invalidated if any of those ownership contracts change.

## Collective Contract

Schedule construction must collectively reject:

- missing or duplicate global source IDs;
- missing or duplicate destination IDs;
- source/destination global ID-set disagreement;
- invalid IDs, integer overflow, allocation failure, or MPI failure.

The forward and reverse paths preserve the caller's local ordering.  All ranks
enter collectives in the same order, including error paths.

## Memory And Scaling

Persistent storage is linear in the locally owned grid count: IDs, owner/position
maps, counts, displacements, and permutation indices.  Per-call storage is two
complex tiles of at most `16 * local_grid_count`.  No `global_grid_count *
orbital_count` array is created.

Schedule construction performs bounded metadata communication once.  Each
Hamiltonian application performs two value-only `MPI_Alltoallv` operations per
tile, plus the existing halo and nonlocal-projector communication inside `hpsi`.

## Verification

Add an MPI test with deliberately permuted and cross-rank source ownership.  For
MPI 1, 2, 4, and 8, require:

- forward redistribution equals a dense global-ID oracle;
- reverse redistribution restores the original values and ordering;
- cached-schedule reuse gives identical results;
- the tiled projected Hamiltonian equals the dense reference;
- duplicate, missing, and out-of-range IDs fail collectively;
- workspace receipts remain proportional to local grid size and tile width.

Keep the route-level contract test that forbids replicated full-cell orbital
storage and requires the production callback to use the cached redistribution.
