# Symmetry-Adapted LCFO Occupied-Subspace Design

## Problem and evidence

The undisplaced Si64 ground state has the expected full atomic symmetry, but
the fragment LCFO truncation does not preserve that symmetry exactly at buffer
boundaries.  The production fixed-center DMN therefore exposes a real closure
failure which the earlier affine proof failed to detect:

- the full atomic catalog contains 1536 affine operations, 32 translations,
  and a 48-operation point cogroup;
- the maximal inversion-containing fixed-center site group has order 12 and
  center `(1/16,1/16,1/16)` for this origin convention;
- the original 384-state LCFO seed gives, for fixed-center operation 4,
  `sigma_min=1.23562e-2`, `sigma_max=9.84658e-1`, and closure residual
  `9.72719e-1`;
- replacing only the 256 unoccupied states by the complete local `s+p`
  projector set improves the residual to `6.65699e-1`, but does not close the
  space.

The remaining leakage is therefore in the 128 occupied LCFO states.  The
previous proof's exactly zero residuals and zero measured workspace are false
evidence and must not be used for publication.

## Decision

Construct a symmetry-adapted occupied subspace once after normal DC
LCFO+EigenExa and before Wannier90.  The construction uses the full-system
atomic symmetry, but stores and transforms spatially distributed fragment
data.  It does not construct or retain the complete full-system orbital
tensor.

For the selected inversion-containing fixed-center group `G`, generate the
orbit of the 128 occupied LCFO orbitals in bounded orbital tiles.  Define the
group-averaged occupied projector

`Pbar = |G|^-1 sum_g U_g P_occ U_g^dagger`.

Select an invariant 128-dimensional subspace from `Pbar`.  Selection is by
complete eigenvalue clusters and group-invariant blocks: a cutoff may not
split a numerically degenerate cluster or an irreducible block.  If no unique
128-dimensional invariant selection exists within the declared tolerances,
reject the checkpoint instead of choosing an MPI-order-dependent basis.

The final 384-dimensional seed is the direct sum of:

1. the symmetry-adapted occupied space of rank 128; and
2. the complete buffer-composed local `s+p` projector space of rank 256.

Orthonormalize the direct sum with the existing generalized algebra.  Reject
rank loss, occupied/projector overlap that prevents rank 384, a non-real
Gamma gauge, or fixed-center closure failure.  No polar-unitary replacement
may hide leakage, no tolerance may be relaxed to accept it, and no
unconstrained Wannier90 or custom-localizer fallback is permitted.

## Symmetry layers

The full affine group and the fixed-center group have different jobs.

- The full atomic catalog defines the physical space-group provenance,
  translation subgroup, point cogroup, origin convention, and fragment
  redistribution rules.  Its operations may have centers outside a fragment.
- The maximal inversion-containing fixed-center group defines the finite
  representation used to adapt the occupied space, write the DMN, constrain
  Wannier90, project operators, and test inversion selection rules.

Both are derived from the full atomic coordinates.  Fragment membership is
used only for distributed storage, buffer support, and returning each MLWF to
the fragment containing its center `R`; it never defines the symmetry group.

## Distributed memory procedure

Spatial-core ranks keep their existing fragment-plus-buffer grid slabs.  For
one group operation and one orbital tile at a time they:

1. pull back the local grid values with the full-system affine operation;
2. exchange only the remote rows required by that pullback;
3. accumulate local contributions to the orbit Gram matrix and to the
   action of `Pbar`;
4. release the transformed tile before advancing.

The eigensolver acts on distributed dense matrices whose dimension scales
with the orbit rank, not the real-space grid.  Production must never retain
`|G|*Nocc*Ngrid`, `Naffine*Nocc*Ngrid`, or a replicated full-system
wavefunction array.  Measured current and peak bytes are reduced across MPI
ranks and recorded; zero measured workspace is invalid when nonidentity
operations are processed.

## Density and boundary receipts

Symmetry adaptation intentionally changes the LCFO projector slightly.  It
must demonstrate that it removes fragment-boundary inconsistency rather than
silently changing the ground state.  Compare the original and adapted
occupied projectors using:

- global occupied-subspace distance and electron-count drift;
- integrated density difference on core-interior grid points;
- integrated density difference in the buffer/boundary shell;
- fixed-center symmetry residual before and after adaptation; and
- minimum selected eigenvalue, maximum rejected eigenvalue, cluster gap, and
  selected invariant-block dimensions.

Interior distortion has the strict tolerance.  Boundary distortion has a
separate, looser diagnostic tolerance because that is where LCFO stitching
error is expected, but remains bounded and published.  These are qualitative
small-cell acceptance checks, not claims of bulk quantitative accuracy.

## Repairing the affine proof

Replace the false-zero proof path with a streamed residual calculation on the
actual selected production seed.  For every required operation it must
measure `||(I-QQ^dagger)U_g Q||`, singular-value bounds of `Q^dagger U_g Q`,
and real allocated workspace.  A fixture that perturbs one operation outside
the seed space must fail.  The proof cannot infer closure merely from atom or
projector permutations.

## Wannier90 and publication

Use the adapted 384-state orthonormal seed as both the band and initial
Wannier gauge for the streamed fixed-center DMN.  Wannier90 remains in
site-symmetry mode and must converge before its iteration limit.  MLWF centers
must close under the fixed-center group, after which orbitals are redistributed
by their centers to owning fragments with buffer support.

V3 records the full affine fingerprint/orders, fixed-center center/order and
fingerprint, before/after closure, density-distortion receipts, invariant
cluster receipt, measured workspace, DMN receipt, Wannier90 convergence, and
center closure.  Any missing, zero where work occurred, nonfinite,
rank-inconsistent, or over-tolerance field rejects restart publication.

## Acceptance sequence

Each implementation task starts with a genuine RED and runs focused MPI
1/2/4/8 verification, specification review, code-quality review, resolution
of every Critical/Important finding, and a clean-first committed-parent
prerequisite overlay build.  The final gate is genuine ideal undisplaced Si64
GS.  Only after it passes may the polarization-primary linear response and
long-pulse HHG run resume; the current remains secondary, and the semilog
polarization spectrum must explicitly classify H2/H4 as peaks or dips and
measure even-harmonic suppression.
