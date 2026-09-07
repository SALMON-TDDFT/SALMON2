# Wannier90 Gamma Symmetry Initialization Design

## Problem

Wannier90 selects `overlap_project_gamma` whenever `gamma_only` is true.  That
routine replaces the complex projection matrix by its real part and performs a
real SVD.  For symmetry-adapted calculations the supplied Wannier
representation can contain non-real translation characters, so discarding the
imaginary part destroys the initial intertwining relation.  The Gamma-specific
routine also omits the `sitesym_symmetrize_u_matrix` call used by the general
complex path.

## Design

Retain the existing real Gamma initialization for ordinary calculations.  When
`lsitesymmetry` is true, select the general complex `overlap_project` path even
at Gamma.  The same condition must be applied in `wannier_lib.F90`, which is
the actual SALMON library entry point, both when choosing `overlap_project`
and when choosing `wann_main`.  That path preserves the complex projection
matrix, projects the initial unitary onto the supplied band/Wannier
representations, and uses the synchronized complex search direction.

No SALMON-side gauge repair is added: the symmetry constraint must hold inside
Wannier90 throughout localization.

## Verification

Add a source-patch regression for the branch condition, rebuild Wannier90, run
the W90 MPI and DMN focused suites, and rerun Si64 through the post-Wannier
generator covariance gate.
