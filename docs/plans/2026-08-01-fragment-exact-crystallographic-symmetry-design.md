# Exact Crystallographic Symmetry for Buffered Fragments

## Purpose

Restore physically correct tensor covariance in the accepted overlapping-Wannier
GS to V3 checkpoint to generalized-eigenvalue Exp coefficient-RT route without
assuming that every material is cubic or that the complete simulation cell is
exactly symmetric.  Symmetry is detected and enforced independently for each
buffered fragment.  The implementation covers nonmagnetic crystallographic
symmetry (32 point groups and operations arising from all 230 space groups).

The normal SALMON route and normal DC LCFO plus EigenExa remain unchanged.  No
removed experimental DG route is restored.

## Physical policy

Only symmetry of the instantaneous atomic structure is admissible.  The
tolerance used while matching atoms, periodic images, and real-space grid
points covers floating-point and coordinate-conversion error only.  It must not
absorb zero-point or thermal displacement.

Consequently, a displaced finite-temperature snapshot may have a smaller
fragment group than its ideal parent crystal, often the identity group C1.  In
that case no nontrivial symmetry projection is applied.  This preserves
physical symmetry breaking, selection-rule relaxation, and dephasing.  Parent
crystal symmetry may be used only to enumerate candidates; it must never be
forced onto the instantaneous structure.  Statistical recovery of parent
symmetry belongs to an ensemble average over independently propagated
snapshots, not to a single-snapshot projection.

## Architecture

### Crystallographic operation source

Reuse the existing optional spglib backend to obtain nonmagnetic space-group
operations for the supplied lattice and atomic structure.  Keep an auditable
operation catalog containing integer fractional-coordinate rotations,
fractional translations, Cartesian rotations, determinant, order, and the
associated crystallographic point-group class.  Tests survey all 32 point
groups and representative nonsymmorphic translations from the 230 space-group
catalog.  SALMON does not maintain a second hand-written symmetry table whose
conventions could diverge from spglib.

An explicit operation-file mode remains available for reproducibility.  Both
sources pass through the same normalization, duplicate rejection, inverse, and
group-closure validation.

### Per-fragment exact subgroup

For each buffered fragment, filter the candidate operations against:

1. atom species and instantaneous positions, including periodic images;
2. the core plus buffer domain and its ownership convention;
3. the discrete real-space grid; and
4. the retained overlapping-Wannier centers and buffered support.

Record separate maximum and RMS residuals for atom, boundary, grid, and center
maps.  Reject operations individually when any residual exceeds its numerical
tolerance.  From the accepted operations select the largest closed subgroup
containing the identity; use a deterministic lexicographic tie break.  Never
repair an invalid operation by moving atoms or increasing the matching
tolerance.  The resulting group may differ between fragments.

Tolerance defaults are scale-aware, based on machine precision, coordinate
roundoff, lattice condition, and grid spacing.  A user may tighten them but
cannot select a thermal-displacement tolerance in this route.  Diagnostics
must make a C1 fallback explicit rather than treating it as a warning or
failure.

### Retained-basis representation and projection

Construct the representation of every retained fragment operation from
buffered real-space overlaps.  Apply polar correction only to remove numerical
non-unitarity, then require metric unitarity, multiplication-table closure, and
inverse consistency.

Use Reynolds averaging within each fragment's exact group:

- scalar matrices: `S` and `H` transform invariantly;
- polar-vector matrices: `X` and `V` transform with the Cartesian rotation;
- improper rotations use the polar-vector convention, including inversion;
- no averaging is performed for C1.

Validate pre- and post-projection covariance.  Projection may remove numerical
noise only; a large pre-projection defect is a hard failure because it signals
an inconsistent basis, operator, or map.  Fragment contributions are assembled
after their local projections, so no global point group is assumed.

### V3 checkpoint and RT

Keep the V3 checkpoint schema.  Fold a deterministic digest of the per-fragment
operation sets, multiplication tables, tolerances, and covariance residuals
into the existing symmetry/observable fingerprint and manifest evidence.
Checkpoint reuse rejects any mismatch in instantaneous geometry or fragment
symmetry data.

Coefficient RT consumes only the projected V3 matrices.  Polarization is the
primary observable and all linear-response and HHG spectra are computed from
field-off-subtracted polarization.  Current remains secondary and is checked
against the centered time derivative of polarization.

## Failure handling

Identity-only symmetry is a valid physical result.  Missing identity,
non-closed accepted operations, non-bijective atom/grid/center maps, excessive
representation defects, inconsistent fingerprints, and non-covariant
post-projection operators are collective fatal errors.  Diagnostics identify
the fragment, operation, residual kind, tolerance, and rejected mapping.

If SALMON is built without spglib, automatic crystallographic discovery is
unavailable and the route must require a validated explicit operation catalog;
it must not silently substitute cubic operations.

## Verification

RED tests precede every implementation unit.  Focused verification includes:

- a catalog survey covering all 32 crystallographic point groups and
  representative symmorphic and nonsymmorphic space-group operations;
- triclinic through cubic lattices, improper rotations, inversion, multiple
  species, periodic images, and fragment-boundary rejection;
- exact ideal structures with only roundoff-scale perturbations;
- physically displaced zero-point/finite-temperature fixtures that reduce to
  the correct subgroup or C1 without projection;
- deterministic per-fragment groups and matrices on 1, 2, 4, and 8 MPI ranks;
- scalar invariance and polar-vector covariance before and after projection;
- unchanged V3 restart identity and forbidden-route contracts;
- genuine ideal Si64 axis equivalence, odd-harmonic dominance, amplitude
  scaling, field-off stability, `J` versus `dP/dt`, and time-step convergence;
- displaced Si fixtures where broken cubic/inversion selection rules are
  reported as physics rather than rejected, while numerical conservation and
  convergence remain mandatory; and
- unchanged normal DC LCFO plus EigenExa regression.

Every task retains RED evidence, focused verification, specification review,
code-quality review, resolution of all Critical and Important findings, and a
clean-first parent-prerequisite overlay build before acceptance.

