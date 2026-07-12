# Scalable SAWF basis-generation contract

This contract is the admission gate for a fixed Wannier sector entering DG-DC.
It applies to gapped, integer-occupation LDA calculations.  `site_symmetry=true`
is a Wannier90 request, not an acceptance criterion.

## Construction routes and symmetry

The small validation cell uses a **monolithic global SAWF** calculation as the
reference.  The hierarchical route classifies representative environments
using the **actual supercell symmetry group**.  For a defective cell the
parent-crystal group is provenance only: a parent-crystal operation absent from
the actual cell is rejected and is never used to restore local symmetry.

The selected bands must be symmetry closed and projections must contain every
member of each complete projection shell.  For every group product the
`D_band` and `D_wann` matrices must obey their representation laws and the
Bloch/Wannier covariance equations within the configured tolerance.  Identity
fallback is forbidden when a nontrivial group was requested.

One SAWF calculation is made for each representative environment orbit.  A
bulk template may be translated/rotated with the validated `D_wann`.  Any core
or buffer intersecting a symmetry-inequivalent defect uses defect-local regeneration.
Independently generated neighbors are joined by buffered
overlap unitary Procrustes alignment before WW, WP, or DG face block assembly.
Rank loss or a gauge residual above tolerance is fatal on the detecting rank,
before a collective reduction.

## Template provenance and cache invalidation

Every template records geometry, pseudopotential, grid, band-window, complete
projection-shell, actual symmetry operations, buffer, schema/code version,
centers, spreads, `D_band`, and `D_wann` fingerprints.  A missing or unequal
fingerprint forces regeneration; stale cache reuse is never silent.  The
recorded gauge unitary and residual form part of the checkpoint.

## Buffer convergence and global reference gate

At least three increasing buffers are required.  The largest two must satisfy
the documented buffer convergence tolerances for centers, occupied projector,
neighbor overlap, WW and WP blocks, and every DG face block.  Localized orbital
convergence alone is insufficient.

On the small cell, one global unitary alignment must bring the hierarchical
basis into agreement with the monolithic reference for occupied projector,
overlap, fixed `H_kin+DG+V_ion`, and every face block.  Only after this
operator-equivalence gate may hierarchical SAWF data initialize DG-DC.

## Namelist controls

All result-changing choices live in `&dc`: `wannier_sawf_generation`
(`monolithic` or `hierarchical`), `wannier_sawf_symmetry_scope` (`actual`),
`wannier_sawf_parent_symmetry_file` (diagnostic provenance only),
`wannier_sawf_buffer_steps` (three increasing integers), gauge, buffer, and
equivalence tolerances, and `wannier_sawf_cache_directory`.
