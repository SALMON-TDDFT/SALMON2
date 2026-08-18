# Distributed Total-Density Buffer Design

## Problem

The stitched overlap-density gate currently receives fragment-local `rho_s` values.  A smooth partition of unity prevents duplicate spatial coverage, but it cannot make independently solved fragment densities equal to the assembled total-system density.  In Si64 the fragment stitch integrates to 254.817110 electrons while `dc%rho_tot_s` integrates to 256.000000000 electrons.

The subsequent ground-state initialization already treats `dc%rho_tot_s` as the physical initial density.  The metric gate must use that same density rather than a different fragment approximation.

## Design

Add an MPI primitive that materializes a row-distributed real scalar field on an arbitrary local list of global physical-grid IDs.  Each rank sends requests for its buffer IDs to the ranks owning the corresponding total-grid slab.  Owners return the scalar values with `MPI_Alltoallv`.  The result is one value per requested local buffer point, in the original request order.

The primitive validates collectively:

- agreed global grid dimensions and tolerance;
- request IDs in `1:N`;
- exactly-once ownership of every global grid ID by the supplied distributed rows;
- finite source and returned values;
- checked default-integer MPI counts, displacements, and byte accounting;
- collective allocation and MPI failures, with cleanup.

Production constructs row IDs for its local `dc%rho_tot_s` slab, calls the primitive once, and uses the returned buffer density in `assemble_dg_stitched_overlap_density_rows`.  It no longer samples fragment `rho_s` for this gate.  No global `N`-element density is replicated.

## Tests

An MPI 1/2/4/8 fixture compares direct redistribution with a dense reference for overlapping and reordered request lists.  Adverse cases cover missing/duplicate ownership, out-of-range IDs, rank-disagreeing metadata, and nonfinite input.  The route checker requires production to use the distributed total-density materialization and forbids fragment-density assignment to `ow_box_density`.

Finally, rerun Si64 MPI 8 / OMP 1 and require the stitched electron count to match 256 before inspecting the subsequent Hermiticity and positive-rank gates.
