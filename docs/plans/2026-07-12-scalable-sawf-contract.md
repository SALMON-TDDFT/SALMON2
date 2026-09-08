# Scalable SAWF basis-generation contract

This contract is the admission gate for a fixed Wannier sector entering DG-DC.
It applies to gapped, integer-occupation LDA calculations.  `site_symmetry=true`
is a Wannier90 request, not an acceptance criterion.

The first DG milestone has two non-negotiable gates: the perfect-crystal
Wannier space must preserve the actual-supercell symmetry through validated
`D_band` and `D_wann`, and DG-DC must produce initial states by solving the
complete DG Hamiltonian including interface, flux, and penalty terms.  These
gates take priority over broader size-scaling studies.

## Construction routes and symmetry

The small validation cell uses a **monolithic global SAWF** calculation as the
reference.  Its symmetry-closed band space and AMN/MMN data come from one
coherent LCFO eigensystem of the same actual supercell.  A conventional
checkpoint is not an accepted SAWF reference.  The hierarchical route classifies representative environments
using the **actual supercell symmetry group**.  The supported non-bulk
environments include defects, material interfaces, free surfaces and vacuum
boundaries, and amorphous regions.  A parent-crystal group is provenance only:
an operation absent from the actual supercell is rejected and is never used to
restore local symmetry.

`wannier_sawf_initial_wavefunction_directory` is reserved for a future
pre-diagonalization seed.  A nonblank value currently fails closed: replacing
LCFO eigenvectors after diagonalization would make the orbitals, eigenvalues,
and AMN/MMN data incoherent.

The selected bands must be symmetry closed and projections must contain every
member of each complete projection shell.  For every group product the
`D_band` and `D_wann` matrices must obey their representation laws and the
Bloch/Wannier covariance equations within the configured tolerance.  Identity
fallback is forbidden when a nontrivial group was requested.

The DC approximation is not repaired by symmetrizing its density or effective
potential after the fact.  `D_band`/`D_wann` unitarity, group closure, and AMN
covariance retain the strict symmetry tolerance.  Hamiltonian covariance is
reported separately and may use the explicit namelist value
`wannier_sawf_hamiltonian_tolerance` to admit the measured DC truncation error.
The accepted residual and configured tolerance are both recorded; this option
must not modify the potential or the orbitals.

For Si, the default valence Wannier sector is the complete atomic `s+p` shell,
four Wannier functions per Si atom.  The `d` shell is not included merely to
enlarge the trial space; `s+p+d` is a separate, explicitly requested convergence
test.  Diamond C uses the same four-per-atom valence-shell default.

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

For the present C64 implementation smoke test, `num_rgrid_buffer=8` is the
accepted code-path value: it exercises symmetry closure, BPW publication, and
all six DG face traces without periodic point duplication.  This smoke-test
choice is not quantitative buffer-convergence evidence and does not relax the
three-buffer gate for later production admission.

On the small cell, one global unitary alignment must bring the hierarchical
basis into agreement with the monolithic reference for occupied projector,
overlap, fixed `H_kin+DG+V_ion`, and every face block.  Only after this
operator-equivalence gate may hierarchical SAWF data initialize DG-DC.

This small-cell gate diagnoses construction and operator consistency, not the
ultimate accuracy of DC-LCFO in the large-system regime.  A small cell can be
dominated by fragment-boundary and buffer fractions.  Existing approximately
4000-atom silicon and water experience, where DC-LCFO wavefunctions reproduced
full-calculation high harmonics through seventh order, is retained as
production-route evidence.  System-size scaling is deferred until the DG
Hamiltonian eigenseed and symmetry-preserving Wannier gates are complete.

SAWF generation and DG-DC admission are deliberately separate phases.  The
generation phase may publish fingerprinted local bases before operator blocks
exist.  It does not thereby authorize their use in DG-DC.  The only
hierarchical handoff is `admit_sawf_hierarchical_basis`, which requires a
matching operator-complete acceptance checkpoint and re-evaluates the
configured tolerance.  A missing, stale, incomplete, or failing checkpoint is
fatal at the handoff; there is no fallback to orbital-only convergence.
The checkpoint carries separate supercell and operator fingerprints.  The
operator fingerprint covers the basis schema, WW/WP definitions, DG face
orientation and penalty convention, fixed-H components, and diagnostic schema;
both fingerprints must match at admission.

## Namelist controls

All result-changing choices live in `&dc`: `wannier_sawf_generation`
(`monolithic` or `hierarchical`), `wannier_sawf_symmetry_scope` (`actual`),
`wannier_sawf_global_reference_source` (currently `lcfo` only),
`wannier_sawf_parent_symmetry_file` (diagnostic provenance only),
`wannier_sawf_buffer_steps` (three increasing integers), gauge, buffer, and
equivalence tolerances, and `wannier_sawf_cache_directory`.

`wannier_sawf_structure_class` declares the intended structure class and may
be `auto`, `crystal`, `defect`, `interface`, `surface`, or `amorphous`.  It
never relaxes exact fingerprint or actual-group checks.  Instead it limits
reuse and records intent: `crystal` rejects unexplained inequivalent
environments, `surface` requires density-measured vacuum occupancy, and
`amorphous` permits reuse only for exact actual-group equivalence. `defect` and
`interface` are non-rejecting policy hints because topology cannot be inferred
reliably from orbit multiplicities alone. No class may relax scientific gates.
