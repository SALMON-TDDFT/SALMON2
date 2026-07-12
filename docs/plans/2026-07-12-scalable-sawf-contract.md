# Scalable SAWF basis-generation contract

This contract is the admission gate for a fixed Wannier sector entering DG-DC.
It applies to gapped, integer-occupation LDA calculations.  `site_symmetry=true`
is a Wannier90 request, not an acceptance criterion.

## Construction routes and symmetry

The small validation cell uses a **monolithic global SAWF** calculation as the
reference.  The hierarchical route classifies representative environments
using the **actual supercell symmetry group**.  The supported non-bulk
environments include defects, material interfaces, free surfaces and vacuum
boundaries, and amorphous regions.  A parent-crystal group is provenance only:
an operation absent from the actual supercell is rejected and is never used to
restore local symmetry.

The selected bands must be symmetry closed and projections must contain every
member of each complete projection shell.  For every group product the
`D_band` and `D_wann` matrices must obey their representation laws and the
Bloch/Wannier covariance equations within the configured tolerance.  Identity
fallback is forbidden when a nontrivial group was requested.

One SAWF calculation is made for each representative environment orbit.  A
bulk template may be translated/rotated with the validated `D_wann`.  Any core
or buffer with a symmetry-inequivalent local-environment fingerprint uses
independent local regeneration.  This includes defect neighborhoods,
interfaces, surfaces/vacuum boundaries, and normally every distinct amorphous
environment.  Independently generated neighbors are joined by buffered
overlap unitary Procrustes alignment before WW, WP, or DG face block assembly.
Rank loss or a gauge residual above tolerance is fatal on the detecting rank,
before a collective reduction.

## Template provenance and cache invalidation

Template reuse is restricted to one exact supercell and its restart lineage;
cross-supercell template import is forbidden.  Every template records the
supercell fingerprint and local-environment fingerprint plus geometry,
pseudopotential, grid, band-window, complete
projection-shell, actual symmetry operations, buffer, schema/code version,
centers, spreads, `D_band`, and `D_wann` fingerprints.  A missing or unequal
fingerprint forces regeneration; stale cache reuse is never silent.  The
recorded gauge unitary and residual form part of the checkpoint.

The supercell fingerprint covers lattice and boundary conditions, all atomic
species and coordinates, pseudopotential content, grid, buffer, band window,
projection shells, XC choice, and generator/schema version.  A local fingerprint
covers the core+buffer species and relative coordinates, local metric,
coordination/distances, and vacuum occupancy.  Environment reuse requires both
an exact local fingerprint match and membership in the same orbit of the actual
supercell symmetry group.  Approximate geometry matching is not an acceptance
criterion.  In an amorphous region this conservative rule normally produces
one independently generated representative per local environment.

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
