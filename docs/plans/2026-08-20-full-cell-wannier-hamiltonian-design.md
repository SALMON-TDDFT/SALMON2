# Full-Cell Wannier Hamiltonian Design

## Problem

The Si64 post-Wannier basis and overlap matrix obey the atomic affine symmetry,
but the stitched Hamiltonian does not.  Direct diagnostics show that the
fragment partition itself is not invariant under the nonsymmorphic point
cogroup: the partition-weight covariance defect is about `7.5e-1`, its
gradient defect is about `1.8`, and the DC density produces a local potential
defect above one.  A Hamiltonian assembled from those fragment windows cannot
serve as a symmetry-preserving final operator.

## Architecture

After the Wannier basis has been constructed, stop using fragment partition
weights and partition gradients to build the final Hamiltonian.  Materialize
the Wannier basis on the complete periodic cell and apply SALMON's established
full-cell Hamiltonian components: periodic finite differences, the updated
full-cell local potential, and the total-system nonlocal pseudopotential.

Wannier orbital `a` is logically owned by the MPI rank that owns its center
`R_a`.  Its spatial values remain distributed over the full-cell real-space
slabs.  Owners schedule bounded orbital tiles (initially 8--16 orbitals); all
ranks apply the full-cell Hamiltonian to their slab of the current tile, and
the projected matrix rows are reduced to the center owners.  No rank gathers
all orbitals on the full grid.

The dense full-cell result is the correctness oracle.  A later change may
drop matrix elements outside a measured center/fragment neighbor range, but
only after comparison with this oracle establishes an error-controlled cutoff.

## Data Flow

1. Finish Wannier90, translation-sector alignment, Gamma sewing, and inverse
   character reconstruction as today.
2. Retain the existing center-to-rank ownership map.
3. Reconstruct the initial full-cell density from the occupied Wannier basis,
   rather than copying the asymmetric fragment DC density.
4. Update the total-system Hartree/XC/local potential from that density.
5. For each center-owner orbital tile:
   - expose only the tile's row-owned coefficients;
   - materialize its values on distributed full-cell spatial slabs;
   - apply the existing SALMON full-cell `hpsi` path using `dc%mg_tot`,
     `dc%system_tot`, `dc%ppg_tot`, and the periodic stencil;
   - form distributed overlaps with all retained Wannier functions;
   - reduce only the resulting owned Hamiltonian rows.
6. Publish row-owned `H`, `S`, density, and operator receipts through the
   existing generalized EigenExa and checkpoint path.

## Memory and Complexity

The persistent dense result remains distributed `O(M^2/P)`.  Full-cell working
storage is `O(G*b)` for tile width `b`, rather than `O(G*M)` on every rank.
The exact dense projection costs `O(G*M^2)` total and is intentionally retained
as the reference.  Future neighbor truncation can reduce storage and apply cost
toward `O(M)` when localization data justify a bounded neighbor count.

## Failure Contracts

- Collectively agree tile width, global grid extent, orbital count, center
  ownership, and basis/operator fingerprints before shape-dependent MPI calls.
- Guard default-integer MPI counts and wide byte receipts before allocation.
- Treat allocation and `hpsi` failures collectively and clean partial tile
  storage on every rank.
- Require finite outputs, Hermiticity, electron count, generalized-eigenpair
  residual, and raw affine covariance before publication.
- Do not repair an order-one defect by group averaging.  Averaging may remove
  roundoff-scale residual only.

## Verification

1. A small MPI fixture compares tiled full-cell projection with a direct dense
   reference on 1, 2, 4, and 8 ranks.
2. The fixture uses nontrivial periodic stencil, local potential, and nonlocal
   projector terms, and verifies center-owner row distribution.
3. A decomposition test changes MPI rank count and tile width while requiring
   identical matrix fingerprints and numerical rows.
4. Existing construction, W90, solver, route, and obsolete-route tests remain
   green.
5. Si64 must reach the one-shot generalized solve with raw `H/S/rho` and
   kinetic/local/nonlocal covariance inside the production tolerance.

## Deferred Sparse Path

Record matrix magnitudes by periodic center distance and fragment-neighbor
shell.  Introduce no cutoff in the reference implementation.  A separate
design will select the smallest shell whose discarded norm, spectra,
transition matrices, and symmetry receipts agree with the full-cell oracle.
