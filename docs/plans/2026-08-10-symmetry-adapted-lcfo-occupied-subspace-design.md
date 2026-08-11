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

The first fixed-center adaptation repaired that finite-group defect, but a
fresh ideal-Si64 run exposed a second, more fundamental error.  The code
discarded the buffer and restricted every seed to its uniquely owned fragment
core before applying the full affine generators.  Measured full-affine
residuals were

- adapted occupied rank 128: total `4.93771`, boundary `4.64122`, interior
  `1.72063`;
- complete `s+p` rank 256: total `8.10362`, boundary `7.82957`, interior
  `3.08775`; and
- their orthonormal rank-384 direct sum: total `10.0621`, boundary `9.32013`,
  interior `3.79228`.

The boundary-dominated failure of both independent blocks disproves the
assumption that a core-truncated complete atomic projector space is already
closed under translations.  Fixed-center adaptation cannot repair a
translation that crosses an artificial fragment face.  The previous proof's
exactly zero residuals and zero measured workspace are false evidence, and a
core-first full-affine proof is physically invalid even when it reports a
finite residual.

## Decision

Construct a smooth full-system distributed seed once after normal DC
LCFO+EigenExa and before Wannier90.  Buffer composition must precede every
full-affine action: multiply each fragment-plus-buffer contribution by the
existing smooth partition of unity, sum equal physical-grid IDs, and retain
one distributed owner for each physical grid point.  Fragment cores may not
be used as hard orbital support during symmetry adaptation or its proof.
This constructs the mathematical full-system functions without constructing
or retaining a replicated full-system orbital tensor.

For the selected inversion-containing fixed-center group `G`, generate the
orbit of the smoothly composed 128 occupied LCFO orbitals in bounded orbital
tiles.  Define the
group-averaged occupied projector

`Pbar = |G|^-1 sum_g U_g P_occ U_g^dagger`.

Select an invariant 128-dimensional subspace from `Pbar`.  Selection is by
complete eigenvalue clusters and group-invariant blocks: a cutoff may not
split a numerically degenerate cluster or an irreducible block.  If no unique
128-dimensional invariant selection exists within the declared tolerances,
reject the checkpoint instead of choosing an MPI-order-dependent basis.

Apply the same smooth composition to all 256 complete buffer `s+p` projectors.
The final 384-dimensional seed is the direct sum of:

1. the symmetry-adapted occupied space of rank 128; and
2. the smoothly composed complete buffer `s+p` projector space of rank 256.

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
- The complete affine generating set acts only after buffer composition and
  proves that the physical rank-384 space, including translations between
  fragment owners, is closed.

Both are derived from the full atomic coordinates.  Fragment membership is
used only for distributed storage, buffer support, and returning each MLWF to
the fragment containing its center `R`; it never defines the symmetry group.

## Cocycle-aware point-cogroup correction

A fresh ideal-Si64 run of the first translation-first implementation disproved
the assumption that the 12-operation fixed-center group completes the affine
adaptation.  The translation average was closed to `1.01211720e-13` and the
subsequent fixed-center average was closed to `8.92736449e-13`, but the final
occupied block still had full-affine residual `3.88036128`.  The reason is
group-theoretic: the fixed-center group contains only 12 of the 48 point
cogroup operations, so translations times that group cover at most 384 of the
1536 affine operations.

The corrected adaptation has three distinct layers:

1. average the smoothly composed occupied projector over the 32-operation
   translation subgroup;
2. average that translation-invariant projector over the 48 affine coset
   representatives; and
3. use the 12-operation inversion-containing fixed-center group only for the
   finite DMN representation, Wannier90 constraint, and inversion receipts.

Coset representatives do not form an exact permutation group: their product
differs from the selected representative by the recorded pure-translation
`global_translation_cocycle`.  Validation must therefore compare a composed
representative action with the point-product representative followed by the
cocycle translation.  It must reject an invalid cocycle, a representative
outside the affine catalog, or a cocycle action inconsistent with the
full-affine product table.  It may not pretend that the 48 representative maps
realize the point-product table exactly.

Because the input projector is already translation invariant, averaging over
one representative from every coset is mathematically identical to averaging
over all 1536 affine operations.  The orbit dimension is therefore
`128*48=6144`, rather than `128*1536`.  The Gram construction should exploit
relative cosets and their translation cocycle so that only 48 distinct
occupied-overlap blocks are evaluated and routed.  Dense EigenExa storage
remains distributed, integer products are overflow-checked before allocation,
and the implementation records separate translation, point-cogroup, and final
full-affine workspace/closure receipts.

Acceptance requires the selected rank-128 projector to pass the existing
strict full-affine generator residual without tolerance relaxation.  A
fixed-center-only closure is diagnostic evidence, not sufficient acceptance.

## Distributed memory procedure

Spatial ranks keep their existing fragment-plus-buffer slabs as input.  A
streamed composition stage routes `(physical_grid_id, weighted orbital tile)`
records to deterministic physical-grid owners.  Contributions from every
overlapping buffer are summed there.  It stores only one orbital tile on the
owned physical grid plus bounded exchange buffers; neither the complete
fragment overlap tensor nor a replicated global grid is retained.

For one group operation and one orbital tile at a time the physical-grid
owners then:

1. pull back the smoothly composed values with the full-system affine
   operation;
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

Only after Wannier90 has produced the global MLWF gauge are orbitals
materialized back onto fragment-plus-buffer slabs.  Each MLWF is assigned by
its center `R` to the fragment containing that center, while neighboring
buffers receive the values required by derivatives and stitched operators.

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

Replace the false-zero and core-first proof paths with a streamed residual
calculation on the smoothly composed physical-grid seed.  Select a
deterministic generating set from
the complete affine product table; invariance under those generators proves
invariance under every generated affine operation, while avoiding a redundant
1536-operation orbital sweep.  For every selected generator it must
measure `||(I-QQ^dagger)U_g Q||`, singular-value bounds of `Q^dagger U_g Q`,
and real allocated workspace.  Full atomic/cocycle provenance still covers
all affine operations.  A fixture that perturbs one generator outside
the seed space must fail.  The proof cannot infer closure merely from atom or
projector permutations, and it cannot use fragment-core truncation as the
support of an orbital.

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
