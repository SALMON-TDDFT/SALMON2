# Wannier90 MLWF Global-LCFO Design

## Decision

Use Wannier90 3.1 as the localization optimizer for the overlapping-Wannier
(OW) ground-state route.  SALMON remains responsible for constructing the
full-system, symmetry-closed LCFO subspace, distributed ownership, physical
validation, fragment assignment, buffer tails, and V3 publication.

The existing normal DC LCFO+EigenExa route and its existing file-oriented
Wannier90/SAWF facilities are unchanged.  The accepted OW route must not
silently fall back to the hand-written occupation-block localizer.

## Why this boundary

The LCFO orbitals already reproduce the self-consistent occupied density and,
at Gamma, have a real gauge.  Localization is therefore a unitary gauge choice
inside the accepted LCFO space, not a reconstruction of the density.  MLWFs
reduce real-space tails and make the later fragment/buffer operator truncation
better conditioned, but minimization alone does not prove crystallographic
covariance.  SALMON must remeasure that covariance after applying the
Wannier90 transform.

Passing full real-space wavefunctions to a global optimizer would recreate the
memory failure.  Wannier90 library mode instead consumes only band-space
matrices and returns unitary matrices, centers, and spreads.  This keeps global
grid values under the distributed SALMON ownership established by the
memory-bounded design.

## Data flow

1. Normal DC-SCF and LCFO+EigenExa produce the authoritative Gamma-real LCFO
   coefficients, occupations, and density receipts.
2. SALMON closes the requested LCFO space under the affine operations derived
   from the full instantaneous atomic configuration.  Fragment geometry does
   not define or modify the symmetry catalog.
3. Distributed point tiles form the reciprocal-neighbour overlap matrices
   `M(m,n,b)` and projection matrix `A(m,w)` required by Wannier90.  Only these
   band-space matrices, eigenvalues, lattice, atoms, and k-neighbour metadata
   are gathered to a designated optimizer rank.
4. The optimizer rank calls `wannier_setup` and `wannier_run` from
   `libwannier.a` with `gamma_only=.true.`.  A minimal generated `.win` controls
   the library and its ordinary `.wout`/checkpoint diagnostics; these files
   are not array spill storage and are never used to move wavefunctions.
5. SALMON validates the returned `U_matrix`, optional disentanglement matrix,
   centers, spreads, and total spread, then broadcasts or block-scatters the
   transform.  The transform is applied to distributed orbital/grid tiles.
6. SALMON recomputes the occupied projector and density, affine cocycle and
   metric covariance, center-orbit closure, Hamiltonian/position/velocity
   covariance, and stationarity.  Wannier90 success is not an acceptance
   receipt.  Any failed receipt collectively rejects V3.
7. Each MLWF is assigned by its periodic center to exactly one fragment core.
   Only support required by configured buffers is redistributed as tails.

## Symmetry policy

The symmetry source is the full-system atomic arrangement and its affine
origin, including centers outside individual fragments.  MLWFs may permute and
mix within a symmetry orbit; individual functions are not required to be
invariant.  The accepted object is the complete covariant Wannier space.

Wannier90's symmetry-adapted inputs may be used when the installed library API
can express the required representation, but they are an optimization
constraint rather than the source of truth.  The post-localization SALMON
receipts remain mandatory for every point group and nonsymmorphic operation.
For a thermally displaced snapshot, symmetry is derived strictly from that
snapshot: no ideal-lattice symmetry is restored and no displacement is added
to the present ideal-Si64 acceptance calculation.

## Memory contract and scalability

Wannier90 3.1 library mode is serial.  It therefore introduces a coordinator
workspace proportional to `nntot*Nband**2 + Nband*Nwann + Nwann**2`, but never
to `Nglobal_grid*Nwann`.  Before gathering any matrix, SALMON must compute the
exact byte requirement with overflow checks, compare it with a configurable
hard limit, and fail collectively before allocation when the limit is
exceeded.  The measured coordinator peak and matrix dimensions become V3
provenance.

For Si64 (`Nband=Nwann=384`, one Gamma point), this matrix-only workspace is
small compared with the stopped replicated wavefunction path.  It is an
appropriate immediate backend, but it is not claimed to provide unbounded
orbital scaling.  Systems beyond the configured dense optimizer limit require
a future distributed Wannier90-compatible optimizer; splitting independent
bundles is not an accepted automatic fallback because it changes the global
spread minimum.

All real-space construction, transform application, symmetry measurement,
and fragment redistribution remain distributed and tiled.  Rank-count changes
may change storage ownership but not the canonical transform (after phase and
permutation canonicalization), centers, coefficients, observables, or V3
payload.

## Determinism and gauge canonicalization

Before applying the transform, SALMON orders Wannier functions by symmetry
orbit and periodic center, fixes each Gamma-real column sign from its largest
canonical core-grid component, and resolves degenerate orbit rotations against
the LCFO projection anchors.  Ties use physical global IDs.  Complex leakage
above tolerance, a non-unitary transform, ambiguous unresolved ties, or
rank-dependent canonical hashes is fatal.

## Build and failure behavior

The OW MLWF route is available only when the bundled Wannier90 library is
linked.  Requesting OW V3 without it is an explicit input/build error.  Normal
DC remains usable without OW.  Library setup/run failure, malformed dimensions,
non-finite output, spread increase beyond its contract, memory-limit excess,
or any post-localization physical failure rejects publication collectively.

## Verification

TDD first compares library MLWF output against deterministic small Gamma
fixtures and records a RED failure while the custom localizer is still called.
MPI 1/2/4/8 tests then require identical canonical transforms and physical
receipts, bounded coordinator memory, no global-wavefunction gather, correct
external symmetry centers, and buffer-only redistribution.  The final clean
overlay enables MPI, ScaLAPACK, EigenExa, spglib, and Wannier90, then reruns the
genuine undisplaced Si64 384/128 GS to V3 before polarization-derived LR and
long-pulse semilog HHG acceptance.
