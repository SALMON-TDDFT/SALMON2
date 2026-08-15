# Spectral-density Wannier seed design

## Goal

Generate symmetry-adapted Wannier90 trial channels without assuming that the
localized states are atom-centred, bond-centred, or interstitial.  The seed
construction must cover covalent bonding/antibonding states, ionic
charge-transfer excitations, molecular crystals, and mixed/interstitial
conduction states.

## Decision

Use the occupied density, the complete empty-state density, and continuous
energy moments of the empty-state density to discover spatial localization
basins.  Production seed generation has no internal empty-state window
boundary.  Equal-count windows remain available only as a diagnostic view.

The densities locate regions; they are not themselves used as orbitals.  For a
basin mask `b(r)` and retained frame `Psi`, construct the small, gauge-covariant
localization operator

```
K_b = Psi^H b(r) Psi.
```

Its numerically separated eigenspaces define trial channels.  This remains
invariant under rotations of degenerate input eigenstates.  Unresolved
eigenspaces stay as blocks and are never split using band number, atom number,
or an arbitrary eigensolver gauge.

## Spectral descriptors

For retained eigenpairs `(epsilon_n, psi_n)` and occupations `f_n`, define

```
rho_occ(r) = sum_n f_n |psi_n(r)|^2
rho_q(r)   = sum_n w_q(epsilon_n) |psi_n(r)|^2
```

For production define the complete empty density and normalized moments

```
rho_empty^(k)(r) = sum_empty ((epsilon_n-E_edge)/E_scale)^k |psi_n(r)|^2,
k = 0, 1, 2.
```

`k=0` contains every empty retained state, while `k=1,2` distinguish the
spatial migration of higher-energy states without a hard spectral cut.  A
tolerance-degenerate block has one energy and therefore one moment weight, so
the descriptors are invariant under its internal unitary gauge.

The spatial feature field contains normalized occupied, unoccupied-window, and
occupied/unoccupied-overlap components.  It therefore represents hole-like,
electron-like, and shared regions without classifying them chemically.

The diagnostic window decomposition is descriptive, never selective, and its
sum is checked against `rho_empty^(0)`.  Basin discovery consumes the occupied
density and all empty moments together.  It may not allocate channel rank or
discard retained states by energy window.

## Basin and symmetry construction

1. Find periodic local maxima of the normalized feature field.
2. Merge maxima whose feature vectors and periodic positions agree within the
   numerical tolerance.
3. Grow deterministic periodic watershed basins from the surviving maxima.
4. Apply every validated affine symmetry operation to each basin.
5. Accept only complete symmetry orbits; symmetry-related basins share one rank
   and one unresolved internal-block structure.
6. Form `K_b` one basin at a time from the full retained frame, not from an
   individual energy window, and retain its separated eigenspaces.
7. Allocate channel ranks orbit-by-orbit until the full retained rank is
   represented.  Failure to span the retained space is a hard rejection.

The resulting target action can be a permutation between basins with a dense
unitary block inside an unresolved local eigenspace.  This block action, rather
than a forced atom-centred monomial action, becomes Wannier90's
`d_matrix_wann`.

## Data flow

The retained LCFO eigenvalues, occupations, and row-distributed physical-grid
values already exist before Wannier90 setup.  New construction is inserted
before the current random-projection setup:

```
retained eigenpairs
  -> equal-count spectral windows
  -> streamed spectral densities
  -> periodic symmetry-closed basins
  -> streamed projected basin operators
  -> trial channel blocks and target symmetry representation
  -> deterministic AMN/projections + d_matrix_wann
  -> Wannier90 site_symmetry minimization
```

Only the densities and one basin operator are resident at once.  No
`grid_points x retained_states x windows` tensor and no full affine-operation
tensor may be allocated.

## Failure contracts

Collectively reject:

- nonfinite or rank-disagreeing eigenvalues, occupations, tolerances, or window
  metadata;
- a spectral boundary that cannot be moved outside a degenerate block;
- nonfinite density accumulation or loss of the spectral partition of unity;
- duplicate/missing spatial ownership;
- incomplete or inconsistent basin symmetry orbits;
- a projected basin operator that is non-Hermitian or loses numerical rank;
- a selected channel catalog that does not span the retained space;
- an AMN/DMN covariance or provenance mismatch.

All extent and byte arithmetic is checked in `int64` before allocation.  Every
rank-dependent failure is reduced before any rank returns from a collective
path.

## Testing

The initial numerical primitive tests cover:

- equal-count windows with a degeneracy crossing a nominal boundary;
- smooth weights that are finite, nonnegative, and sum to one;
- invariance under permutation and unitary rotation of degenerate states;
- separated occupied/electron densities for an ionic charge-transfer model;
- bonding/shared and interstitial maxima for a covalent model;
- rank-disagreeing metadata, nonfinite inputs, unsafe extents, and incomplete
  spectral coverage.

Integration tests then verify MPI 1/2/4/8 equality, complete symmetry orbits,
AMN/DMN covariance, Wannier90 convergence, and Si64 point-cogroup closure.  The
existing post-Wannier centre canonicalizer remains a validator during the
transition and is removed only after the native symmetry-adapted path passes.
