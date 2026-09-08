# Memory-Bounded Global Wannier Symmetry Design

## Goal

Keep the exact full-system crystallographic symmetry used to construct the
Gamma-point global LCFO Wannier space, while removing replicated full-orbital
and full-representation storage before fragment redistribution.  The operation
runs once before RT, so peak memory is the primary optimization target; disk
spill is forbidden and additional arithmetic is acceptable.

## Observed scaling problem

The ideal Si64 run reaches the correct 384-state LCFO target and factors its
1536 affine operations into 32 pure translations and 48 point-cogroup
representatives.  However, `measure_dg_rank_fixed_symmetry_residuals` currently
replicates on every rank:

- the complete local orbital rows in `orthonormal_basis`, `owner_basis`,
  `image`, and `residual`;
- every `384 x 384` representation for all 48 representatives; and
- temporary dense overlap and metric matrices.

It also broadcasts every owner's complete orbital block for every operation.
The stopped eight-rank Si64 run measured about 1.99 GiB RSS per rank, or about
15.9 GiB in aggregate, before V3 publication.  This ownership model is not an
acceptable large-system procedure.

## Physical invariants

The memory change must not weaken the physics:

1. Symmetry operations come only from the full instantaneous atomic
   configuration, not from independently symmetrized fragments.
2. Pure translations are factored from point representatives while fractional
   translations and multiplication cocycles remain exact.
3. The retained LCFO space contains the authoritative occupation spectrum and
   all requested columns.  At Gamma, non-negligible imaginary LCFO coefficients
   remain a hard error.
4. Individual Wannier functions need not be invariant and no common fixed
   center is required.  The complete center orbit, metric representation,
   density, Hamiltonian, and observables must be covariant.
5. After construction, a Wannier is assigned by its periodic center to one core
   fragment.  Only its required buffer tails may be replicated.
6. V3 publication remains impossible unless every GS acceptance receipt and
   provenance fingerprint passes.

## Distributed ownership

Use a logical two-dimensional decomposition:

- the existing spatial dimension owns unique core grid points and fragment
  buffers;
- a new orbital dimension owns contiguous or block-cyclic ranges of global
  LCFO/Wannier indices.

No rank may allocate the full `Norb x Nlocal` basis merely to apply symmetry.
Each rank holds only its orbital rows on its spatial points.  Small global
integer metadata, the 48-operation catalog, and scalar residual receipts may be
replicated.

Dense `Norb x Norb` metric, overlap, and representation blocks use the existing
ScaLAPACK-compatible block-cyclic ownership.  A serial dense copy is permitted
only in focused fixtures below an explicit small-size threshold; production
Si64 and larger cases must use the distributed path.

## Streaming symmetry action

Process one point-cogroup representative at a time:

1. Determine remote core values required by the representative's exact point
   permutation.
2. Exchange only those values with `MPI_Alltoallv`; do not broadcast complete
   owner bases.
3. Form the distributed overlap block and reduce it into the distributed
   representation.
4. Accumulate total, boundary, and interior residual norms in point blocks.
   Do not allocate simultaneous full `image` and `residual` arrays.
5. Verify metric unitarity and the required cocycle products for the current
   block.  Update the affine-cocycle fingerprint and maximum receipt values.
6. Release operation-local buffers before advancing to the next representative.

Pure translations are applied through their integer point permutation and
Gamma phase rule.  They are not expanded into 32 separately stored dense
representations.  Product validation uses the point representative plus the
stored translation cocycle.

## Global localization and redistribution

The global localization transform is distributed by Wannier index.  Center and
spread reductions produce only `O(Norb)` replicated metadata.  Orbital values
remain distributed until every periodic center has a deterministic core
fragment owner.  Redistribution then sends each core orbital block to that
owner and sends tails only to fragments whose configured buffer contains the
corresponding support.

The redistribution fingerprint includes the orbital owner map, physical core
IDs, tail IDs, generation counters, and communicator layout.  Rank-count
changes may alter storage ownership but must not alter physical coefficients,
occupations, centers, operators, or checkpoint payload after canonical
serialization.

## Memory contract

Production code records peak allocated bytes for the symmetry/localization
workspace.  The acceptance contract requires:

- no production allocation proportional to `Norb * Nglobal_grid` on one rank;
- no replicated allocation proportional to `Nsym * Norb**2`;
- operation-local workspace released between representatives;
- per-rank large-array memory decreasing when the same fixture is run on
  1/2/4/8 ranks, within small replicated metadata overhead; and
- a documented hard failure when the requested process grid cannot distribute
  a required dense block safely.

RSS remains supporting evidence because libraries allocate outside the module,
but deterministic module allocation counters are the primary regression gate.

## Failure handling

Fail collectively before V3 publication for an invalid process grid, missing
point-map target, inconsistent send/receive counts, non-finite distributed
block, singular metric, Gamma-real violation, failed metric unitarity, failed
cocycle closure, center-orbit failure, uncovered buffer tail, or mismatched
canonical fingerprint.  All ranks must take the same rejection path.

## Verification

TDD fixtures cover distributed ownership, exact Alltoallv point permutation,
streamed versus dense reference residuals, cocycle products, allocation bounds,
and canonical rank-count identity on MPI 1/2/4/8.  Production verification then
uses a clean-first MPI+ScaLAPACK+EigenExa+spglib overlay and the genuine ideal
Si64 384/128 route.  Only after the memory gate, GS receipts, inversion, and V3
publication pass may polarization-derived LR and long-pulse HHG be run.
