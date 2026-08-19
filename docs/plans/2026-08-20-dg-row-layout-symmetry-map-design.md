# DG Row-Layout Symmetry-Map Reindexing Design

## Problem

The post-Wannier DG path changes `ow_core_ids` from the row layout used to build
`global_symmetry_map` and `fixed_center_symmetry_map` to a physical-grid-ID
layout.  The current call passes those physical IDs to
`exchange_dg_point_permuted_orbital_rows`, although that routine expects global
row indices.  The two symmetry maps are then retained in their old row-index
coordinate system.  Values can remain accidentally self-consistent with the
stale maps, while finite differences use the new physical layout; this explains
the observed order-one gradient covariance defect.

## Chosen Approach

Keep the desired post-Wannier physical-ID row layout, but perform the transition
explicitly:

1. Materialize `global_closed_core` on `initial_core_ids` using the existing
   physical-ID keyed `materialize_ow_distributed_core_to_buffer` routine.
2. Reindex every target row in `global_symmetry_map` and
   `fixed_center_symmetry_map` from the old distributed row layout to the new
   distributed row layout by matching physical IDs.
3. Only after both operations succeed, assign `ow_core_ids=initial_core_ids`.

This is preferred over retaining the old layout throughout the remaining
pipeline or rebuilding the crystallographic catalog.  It is local to the actual
coordinate-system transition and preserves all established operation ordering,
catalog fingerprints, and downstream algorithms.

## Reindexing Primitive

Add a distributed helper in `dg_overlapping_wannier_construction` with inputs:

- old local physical row IDs;
- new local physical row IDs;
- one or more old-layout target-row maps.

The helper collectively validates rank-agreed operation count, exactly-once old
and new physical-ID ownership, equal global row counts, target-row ranges, and
checked allocation/MPI extents.  It constructs two temporary global lookup
tables: old global row to physical ID, and physical ID to new global row.  For
each new local source row it finds the corresponding old source row, reads its
old target row, and converts that target through physical ID into the new global
row index.  Output remains row-owned and no dense orbital matrix is introduced.

All allocation and validation failures are reduced collectively before return.
The production caller keeps the old values, IDs, and maps unchanged until every
new object has been created successfully.

## Testing

Use test-driven development.  First add an MPI fixture with deliberately
different old and new row orderings on every rank.  Define a nontrivial physical
permutation and prove that applying the reindexed map in the new layout produces
the same physical action as the original map in the old layout.  Run it on MPI
1, 2, 4, and 8 ranks and verify the test fails before the helper exists.

Add adverse checks for duplicate/missing physical ownership and invalid target
rows.  Update the route test to require physical-ID materialization and map
reindexing before `ow_core_ids` changes, and to forbid the erroneous direct
exchange call at this transition.

After implementation, run the focused construction MPI suite, route and memory
lifetime checks, production build, and then the Si64 run with MPI rank count 8
and `OMP_NUM_THREADS=1`.  The acceptance signal is that the post-reorder
finite-difference/map commutator and gradient covariance defect become small;
the downstream physical validation gates remain unchanged.

## Scope Control

This change does not rebuild symmetry catalogs, alter the Wannier gauge,
average a defective Hamiltonian, change MPI rank count, or add a new physical
approximation.  It only repairs the row-coordinate transition already intended
by the production path.
