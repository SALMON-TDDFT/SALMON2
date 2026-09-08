# Translation-Covariant Post-MLWF Gauge Design

## Context

The global LCFO space is closed under the 1536-element affine symmetry group
to about `1e-12`, and the fixed-center 12-operation DMN now reaches Wannier90.
Wannier90 nevertheless breaks the 32-element atomic translation subgroup:
the resulting center set has a worst nearest translated-center mismatch of
about `0.286` in supercell fractional coordinates.  This is not a tolerance
problem.  Only the fixed-center subgroup is constrained during localization.

The ground-state basis must retain the atomic translations.  RT states must
not be forced to retain them, because a spatially varying external field may
break translation symmetry.  Memory, rather than one-time preprocessing
runtime, is the primary large-system constraint.

## Decision

Use Wannier90 to obtain the initial MLWF gauge, then apply a bounded-memory
translation-orbit gauge correction.  Exact translation covariance takes
priority over the unconstrained Wannier90 minimum.  The correction is the
smallest unitary change within each translation orbit, so it preserves as much
of the MLWF localization as the symmetry constraint permits.

Do not send all 1536 supercell affine operations to Wannier90.  Do not build a
primitive-cell-only route: defects, interfaces, nonprimitive supercells, and
arbitrary fragment layouts must use the same procedure.

## Mathematical construction

Let `T` be the 32-element atomic translation subgroup and let `W` be the
post-Wannier90 basis.  Use the already validated point maps to stream the
overlap action

`B_t = <W | P_t W>`

one translation at a time.  Recover translation orbits from these overlaps
without storing `|T| * Nwann^2` matrices.  Choose a deterministic representative
per orbit from the canonical pre-MLWF ownership and center metadata.

For an orbit representative block `W_0`, generate the target translated block
`P_t W_0`.  Align the corresponding Wannier90 block by the unitary polar factor
of its cross overlap.  Assemble a block-sparse correction `Q_T`.  Apply a final
metric polar correction within each complete orbit so that the corrected basis
is orthonormal and obeys the translation multiplication table to the requested
tolerance.

The correction must not mix the occupied-derived and projection-localizer
rank contracts in a way that changes the retained rank.  It acts on the final
384-dimensional basis only and is applied identically to values, gradients,
operator transformations, and provenance fingerprints.

## Point-cogroup and cocycle

After translation correction, validate the 48 point-cogroup representatives
with the existing translation cocycle.  Products are checked as representative
plus translation, rather than by materializing all 1536 representation
matrices.  This proves full affine closure while retaining `O(Nwann^2/P)`
distributed storage plus one streamed dense operation on the coordinator.

The fixed-center 12-operation DMN remains responsible for Wannier90 site
symmetry.  The post-MLWF translation correction is responsible for the symmetry
that the compact DMN deliberately omits.

## Data flow

1. Build the buffer-supported global LCFO basis and prove full affine closure.
2. Publish the fixed-center DMN and run Wannier90.
3. Apply Wannier90's unitary transform to core and buffer values and gradients.
4. Discover complete translation orbits deterministically.
5. Stream translation images and construct orbit-local polar alignments.
6. Apply the distributed block-sparse gauge correction.
7. Recompute centers and redistribute each orbital to the fragment containing
   its corrected center; the symmetry center may lie outside that fragment.
8. Validate translation products, point-cogroup cocycles, density invariance,
   localization cost, and full affine closure.
9. Write the V3 checkpoint and use the fixed corrected basis for Exp RT.

## Memory model

Never allocate `Ntranslation * Nwann^2`, `Naffine * Nwann^2`, or transformed
real-space copies for every group operation.  Row-owned overlaps remain
distributed.  The root may gather one `Nwann^2` operation at a time.  Orbit
corrections are stored as small dense blocks plus `O(Nwann)` integer metadata.
All allocation and MPI-count arithmetic is overflow checked and reported in
peak-workspace diagnostics.

## Acceptance gates

- Translation orbit decomposition is complete and rank independent on MPI
  1/2/4/8.
- Translation identity, unitarity, and product residuals satisfy the configured
  symmetry tolerance.
- The 48-representative cocycle proof implies the full 1536 affine closure.
- Corrected center orbits close periodically; no tolerance relaxation is used
  to hide a broken orbit.
- Electron count and density are invariant under the unitary correction.
- Rank remains exactly 384 for the Si64 acceptance case.
- The spread increase is measured and bounded by an explicit acceptance
  threshold; both pre- and post-correction values are reported.
- The procedure works for arbitrary fragment shapes and does not infer
  symmetry independently per fragment.
- Normal DC LCFO+EigenExa and Exp-only RT routes remain isolated.

## Verification and review

Each implementation task starts with a genuine RED.  Run its focused MPI
1/2/4/8 verification, specification review, and code-quality review.  Resolve
all Critical and Important findings before committing.  Every task also uses a
clean committed-parent, full-feature overlay build with MPI, OpenMP, EigenExa,
ScaLAPACK, Spglib, and Wannier90.

The genuine ideal Si64 run must pass fixed-center DMN, Wannier90, translation
correction, full-affine center/representation gates, and V3 checkpoint before
linear-response or HHG work resumes.
