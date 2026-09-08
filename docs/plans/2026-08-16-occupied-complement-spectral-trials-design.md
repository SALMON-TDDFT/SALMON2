# Occupied/Complement Spectral Trial Frame Design

## Goal

Construct a complete, symmetry-compatible Wannier90 trial frame without
selecting mutually redundant spectral-basin channels from the full retained
space.  Preserve the already adapted occupied frame and apply spectral-basin
localization only to its orthogonal complement.

## Context

For the Si64 retained space, `ntarget=384` and `nstate=128`.  The current
spectral-basin route constructs occupied and empty density descriptors, but it
projects every basin operator into all 384 retained directions and selects a
total of 384 basin channels independently.  The propagated candidates then
have Gram defect 1.0: the basin channels duplicate directions instead of
forming a complete frame.

The downstream density route does not identify Wannier functions themselves
as occupied orbitals.  It solves for `nstate` physical coefficient vectors in
the final retained basis and applies the material-dependent occupation vector
to those vectors.  Occupied/complement separation is therefore a stable trial
frame construction rule, not a hard physical occupation label on the final
Wannier functions.

## Chosen Architecture

1. Derive dimensions from production metadata:
   `noccupied = nstate`, `nempty = ntarget - nstate`.  Do not use Si-specific
   constants or MPI-rank-dependent counts.
2. Keep the first `noccupied` columns of the already symmetry-adapted retained
   frame as the occupied trial block.  Do not copy or re-diagonalize them
   through every basin.
3. Construct the localization complement from the existing complete-s+p
   projector seeds after their joint orthonormalization against the adapted
   occupied LCFO block.  DC-LCFO diagonalizes `ntarget` states internally, but
   production deliberately requests real-space buffered contributions for only
   `nstate` states; it must not reintroduce all empty LCFO buffer wavefunctions.
   The API validates the occupied/complement Gram and cross-Gram defects rather
   than assuming that a trailing coordinate range proves the physical origin.
4. Project each prepared spectral-basin operator directly into the `nempty`
   complement.  Store only one `nempty x nempty` operator and one eigensystem
   at a time.
5. Select and propagate exactly `nempty` complement channels across basin
   orbits.  The propagated empty block must be orthonormal and orthogonal to
   the preserved occupied block before concatenation.
6. Assemble one distributed `ntarget`-column trial frame
   `[occupied | localized complement]`.
7. Build one AMN/MMN/DMN input set and invoke Wannier90 exactly once over the
   complete retained space.  No Wannier90 call is made per occupation block,
   basin, or basin orbit.

## Memory and Computation

Let `N` be the real-space grid size, `M=ntarget`, `O=nstate`,
`E=M-O`, `B` the number of basin orbits, and `P` the MPI rank count.

The basin projection and diagonalization cost changes from

`O(B*N*M^2/P + B*M^3)`

to

`O(B*N*E^2/P + B*E^3)`.

For Si64, `M=384` and `E=256`, so the projection coefficient is reduced to
about 44 percent and the dense basin eigensolve coefficient to about 30
percent.  The final one-shot Wannier90 work remains approximately
`O(N*M^2/P + I*M^3)` and therefore does not acquire ideal asymptotic weak
scaling.

Per-rank owned memory remains

`O(N*M/P + E^2 + M^2/P)`

plus the existing bounded Wannier90 coordinator workspace.  The implementation
must not allocate all basin operators, separate occupied and empty real-space
frames, or a replicated `N x M` frame.

## Invariants and Failure Handling

- `0 <= noccupied <= ntarget` and `nempty=ntarget-noccupied` on every rank.
- Replicated dimensions and tolerance agree collectively before any
  shape-dependent collective.
- The occupied frame is orthonormal and symmetry closed within tolerance.
- The complete-s+p-derived complement frame is orthonormal, complete, and
  cross-orthogonal to the occupied LCFO frame.
- Basin rank selection returns exactly `nempty`, not `ntarget`, channels.
- The concatenated trial frame has Gram defect within tolerance.
- All extent and workspace arithmetic is checked before allocation.
- Rank-local allocation or validation failure is converted into a collective
  rejection before later MPI calls.
- `nempty=0`, `noccupied=0`, and trivial translation groups are valid boundary
  cases.

## Test Strategy

1. Add a focused failing fixture in which full-space independent basin
   selection duplicates an occupied direction, while complement-only
   selection produces a complete orthonormal frame.
2. Verify a material-dependent non-half split; no test may derive dimensions
   from MPI rank count.
3. Verify empty-only, occupied-only, duplicate-channel, cross-Gram, extent,
   rank-disagreement, and allocation failure paths.
4. Run the construction and EigenExa MPI fixtures on 1/2/4/8 ranks and compare
   decomposition-independent receipts.
5. Run route and obsolete-route checks, overlay build, and `git diff --check`.
6. Re-run Si64 through the previous failure point and require the concatenated
   trial Gram gate to pass before starting the single Wannier90 call.

## Deferred Scaling Work

This change removes unnecessary full-space basin work but does not convert the
overall dense Wannier route into a weak-scaling algorithm.  Sparse/local
orbital-pair overlap, block-local Wannier optimization, and hierarchical
collectives remain separate future work.
