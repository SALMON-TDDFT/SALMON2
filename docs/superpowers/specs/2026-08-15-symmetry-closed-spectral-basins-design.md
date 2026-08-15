# Symmetry-Closed Periodic Spectral Basins

## Problem

The periodic spectral-basin watershed resolves tolerance-equal neighbours by
choosing the smaller global grid ID.  That naming rule is deterministic but is
not equivariant under translations or point operations.  Symmetric flat
regions can therefore be split into a partition on which a supplied symmetry
generator does not induce a basin permutation.  The Si64 production run
reaches this condition after all retained-space symmetry closure checks have
passed at approximately `1e-12`.

The generator maps are row-owned.  The solution must not replicate an
`Npoint x Ngenerator` action on every MPI rank.

## Design

Keep the existing six-neighbour watershed as the preliminary fine partition.
Then compute the finest coarsening of that partition that is invariant under
the supplied generator permutations.

Represent the partition by a replicated union-find over `Npoint` points.  Its
initial equivalence classes are the preliminary watershed basins.  Stream one
generator map at a time into a temporary `Npoint` integer vector.  For every
currently equivalent point pair, require their generator images to remain
equivalent.  Union image classes when this condition is violated.  Repeat over
all generators until a complete pass makes no union.  Because the generators
are finite permutations generating the supplied symmetry group, the fixed
point is the smallest generator-invariant equivalence relation containing the
original watershed partition.  The induced action on the final classes is a
permutation.

Global point IDs may be used to choose union-find representatives and final
label order.  They affect only canonical naming, not which points are merged.

## Distribution and Memory

The public input remains `generator_maps(nlocal, ngenerator)` paired with
`row_ids(nlocal)`.  For each generator, ranks place their owned rows in a zeroed
temporary global vector and combine it collectively.  Exactly-once row
ownership and map range are checked before indexing.

Additional persistent workspace is bounded by a small number of `Npoint`
integer vectors.  No `Npoint x Ngenerator` array is allocated.  Workspace byte
receipts are checked before allocation and reduced with `MPI_MAX`.

## Failure Handling

All allocation, extent, map, convergence, and MPI failures are reduced
collectively before return.  Closure iteration is bounded by `Npoint - 1`
successful class merges; exceeding the bound rejects collectively.  The final
basin action is independently checked for a complete one-to-one permutation.

## Verification

1. Add a flat periodic-grid RED with a nontrivial translation for which the old
   global-ID tie break creates a non-closed partition.
2. Require the new result to be the minimum symmetry-closed coarsening and to
   have a valid generator permutation.
3. Keep the row-distributed-map fixture and compare basin fingerprints across
   MPI 1, 2, 4, and 8.
4. Retain out-of-range and ownership rejection tests.
5. Run construction MPI, EigenExa MPI, route, full build, and the Si64
   production case through the spectral-basin stage.

