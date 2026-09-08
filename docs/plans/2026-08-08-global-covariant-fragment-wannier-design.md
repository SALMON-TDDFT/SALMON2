# Global-Covariant Fragment-Local Wannier Design

## Purpose

Preserve every exact instantaneous full-system symmetry in the retained
overlapping-Wannier ground-state route without requiring a symmetry center to
lie inside a fragment and without gathering full-system real-space
wavefunctions onto one rank.

The accepted output remains the buffer-supported overlapping-Wannier GS, V3
checkpoint, and generalized-eigenvalue Exp coefficient RT.  Conventional
SALMON and normal DC LCFO plus EigenExa remain unchanged.

## Physical definition

Symmetry belongs to the instantaneous full atomic system, not to a fragment.
For a full-system affine operation `g = (R_g, tau_g)`, define one global
Wannier frame by

```text
psi_{g n}(r) = sum_m D_mn(g) psi_m(g^-1 r).
```

`D(g)`, the Wannier-center orbit, periodic lattice translation, and associated
phase convention are fixed once for the full system.  A fragment stores only
the restriction of these globally defined functions to its core and buffer.
The symmetry center may lie outside that fragment.

An operation is not required to map one fragment box to one fragment box.  It
may split the points of one source fragment among several target fragments.
Each grid point is therefore mapped by its global periodic grid ID and
resolved to its actual owner rank.

## Architecture

### Full-system symmetry catalog

Construct the crystallographic catalog only from the full instantaneous
lattice, atomic coordinates, and species.  Retain operations commensurate with
the global real-space grid.  Do not reject an operation because it fails to
preserve a fragment boundary or buffer shape.

At finite ionic displacement, use the exact instantaneous group found at the
configured tolerance.  Never restore a parent operation absent from the
instantaneous structure.

### Global gauge and phase convention

Choose deterministic representative Wannier centers and deterministic
representative orbitals.  For every retained operation record:

- the affine global-grid point map;
- the center-orbit map;
- the orbital mixing or permutation matrix `D(g)`;
- the periodic lattice wrap associated with each mapped center; and
- its lattice-translation phase.

The representation must satisfy metric unitarity and the complete affine group
product table.  Gauge fitting is not performed independently by fragment.

### Fragment-local storage and communication

Each rank retains the Wannier values and gradients needed on its local core and
buffer only.  When an operation maps those points to several owners, the rank
fetches the corresponding pieces owner by owner.  No rank materializes the
full-system real-space grid for every Wannier.

Only bounded point batches and dense Wannier-space matrices are collective.
Memory scales with local buffered grid storage plus `O(N_W^2 |G|)` small-matrix
evidence, not with `N_grid N_KS |G|` on one rank.

### Construction and localization

Build the candidate and retained subspaces using the global pointwise action.
The retained space must contain the occupied subspace and must be closed under
every retained `D(g)` before localization.  Insufficient target rank is fatal;
the implementation must not silently drop an operation or enlarge the basis.

Localization uses transformations commuting with the global representation.
It may evaluate and update fragment-local restrictions, but every update is
tied across the complete global symmetry orbit.  Publication requires both
spread convergence and full-system subspace closure.

### Operator publication

Measure the symmetry action from the actual distributed localized Wannier
frame.  Validate, before any projection:

- subspace leakage;
- metric unitarity;
- affine group closure;
- scalar covariance of `S` and `H`; and
- vector covariance of polarization/position and velocity.

Projection may remove roundoff-scale residuals only.  A physical-scale defect
must reject V3 publication.  In particular, ideal centrosymmetric Si must
publish inversion evidence before its checkpoint can be used for HHG.

## Failure handling

Collective validation failures are rank-consistent and fail closed.  The route
rejects publication when the global operation is grid-incommensurate, a mapped
point has no owner, the target space is not symmetry closed, representation
closure fails, localization stalls, or pre-projection covariance is larger
than the accepted tolerance.

No failure may fall back to a removed WPW, Fragment, Nodal, mixed-z, full-H
seed, or adaptive-DG-basis route.

## Verification

The implementation uses RED-first MPI fixtures covering:

- an operation whose center is outside a fragment;
- one source core split across multiple target ranks;
- periodic wraps and phase/group-product consistency;
- arbitrary fragment numbering and decompositions on 1, 2, 4, and 8 ranks;
- occupied-subspace inclusion and fixed-rank symmetry closure;
- rejection of an insufficient target rank;
- localization that preserves the global representation; and
- ideal Si64 inversion, V3 publication, polarization-derived linear response,
  and a semi-log HHG spectrum with odd peaks and classified even peak/dip
  morphology.

Every implementation task includes focused verification, specification review,
code-quality review, resolution of all Critical/Important findings, and final
clean-first parent-prerequisite overlay verification with MPI, ScaLAPACK,
EigenExa, and spglib enabled and Wannier90 disabled.
