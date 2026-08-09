# Global LCFO Orbital-Parallel Wannier Design

## 1. Purpose

Generate 384 localized real Wannier orbitals from one coherent full-system DC-LCFO eigenspace, then redistribute each localized orbital to the fragment containing its periodic center.  Ground-state density and the occupied one-particle projector must be preserved independently of fragment boundaries and Wannier gauge.

This design replaces fragment-local occupied closure, fragment-local complement generation, and dense enumeration of every supercell space-group operation.

## 2. Physical invariants

### 2.1 Occupied projector, not density alone

The converged Gamma-point DC state has real orbitals and a symmetry-compatible density.  The density diagonal alone does not determine orbital signs or the occupied subspace.  The authoritative ground-state invariant is

\[
P_{\mathrm{occ}}(\mathbf r,\mathbf r')=
\sum_{n=1}^{N_{\mathrm{occ}}}
\psi_n(\mathbf r)\psi_n(\mathbf r'),
\qquad N_{\mathrm{occ}}=128\ \text{for Si64}.
\]

The Wannier construction must preserve the LCFO metric, the occupied projector, electron count, density, and occupied energy within explicit numerical tolerances.

### 2.2 Real Gamma gauge

At Gamma without spin-orbit coupling, the retained LCFO Hamiltonian and orbitals are real.  Gauge freedom is therefore a real orthogonal transformation, including signs and rotations inside degenerate subspaces.  Complex phases are rejected on this route rather than silently discarded.

### 2.3 Occupation-preserving blocks

For an insulating zero-temperature case, the default localization transformation is block diagonal,

\[
U=U_{\mathrm{occ}}\oplus U_{\mathrm{unocc}},
\]

with dimensions 128 and 256 for Si64.  This preserves the occupied projector exactly.  Occupied-unoccupied disentanglement is outside this route; adding it would require a separately specified frozen/outer-window projector-preservation algorithm.

For a general material with fractional occupations, the invariant is the one-particle density matrix `P = C f C^T`, not an integer occupied rank.  Partition the retained LCFO states into blocks with equal occupation eigenvalues within a tolerance fixed before localization.  Gauge rotations are allowed only inside an equal-occupation block.  Equivalently, the localization transform must commute with `f`.  States with unequal fractional occupations may not be mixed.  Spin-polarized real-Gamma support requires the same construction independently for each spin channel; spin-orbit and non-Gamma complex gauges remain outside the retained route.

## 3. Input space

DC-LCFO constructs the full-system Hamiltonian in localized fragment basis functions \(\phi_\mu\).  EigenExa returns at least 384 real coefficient columns,

\[
\psi_n(\mathbf r)=\sum_\mu \phi_\mu(\mathbf r)C_{\mu n},
\qquad n=1,\ldots,384.
\]

Normal DC LCFO+EigenExa output remains unchanged.  The overlapping-Wannier route requests additional coefficient columns in memory when `nstate_tot < 384`.  The request is rejected if the LCFO basis rank is below 384.

Every retained orbital comes from this same coefficient matrix.  Atomic s+p projectors may define a localization objective, labels, and initial gauge, but they never append basis vectors.

## 4. Data distributions

### 4.1 Fragment distribution before localization

LCFO basis functions, basis transforms, and Hamiltonian blocks remain distributed by fragment.  After EigenExa extraction, each fragment root retains its LCFO-basis row block and the requested retained coefficient columns.  The design does not require a new persistent two-dimensional distribution of `C`; bounded column batches may be broadcast within the fragment communicator and exchanged only as coefficient data.

No rank gathers all fragment basis functions or all 384 full-system real-space orbitals.

### 4.2 Orbital ownership during materialization

The 384 Wannier indices are balanced-block distributed over MPI ranks, with owner counts differing by at most one.  With eight ranks, each rank owns 48 materialized full-system orbitals.  MPI ranks therefore have two roles: fragment ranks perform spatially local basis contractions, and orbital-owner ranks assemble bounded batches of complete orbitals.  This is an orbital/spatial transpose, not a migration of the LCFO basis.  OpenMP parallelizes grid points and local basis contractions.

LCFO basis arrays never leave their fragment ranks.  For one bounded orbital batch, distribute only the small transformed coefficient slice.  Every fragment rank contracts its resident LCFO basis on its unique core for that batch.  An `MPI_Alltoallv` transpose sends the disjoint core slabs to the ranks owning those Wannier indices, which assemble their assigned full-system orbitals.  Counts and displacements are checked against the MPI integer range before allocation or communication.  Communication is bounded by `orbital_batch_size`; it does not all-gather all orbitals or basis functions.  Across all batches, moving each materialized orbital once is unavoidable, but the peak memory remains bounded.  Buffer values are not included in this first transpose because unique cores already cover the complete periodic grid.

For Si64 on a 32-cubed grid, 48 real full-system orbitals occupy about 12 MiB per rank.  Temporary storage must remain proportional to `orbital_batch_size * global_grid_count`, not `384 * global_grid_count` on every rank.

### 4.3 Spatial redistribution after localization

For each localized orbital, compute its three periodic first moments and fractional center.  Determine the owning fragment from the unique periodic core partition.  The partition uses periodic half-open cells; an exactly face-, edge-, or corner-centered orbital is assigned by the smallest global core-grid identifier and then fragment identifier.  Geometric center and computational owner are separate metadata; a center may initially be outside the rank that constructed the orbital.

Perform an MPI all-to-all transpose from orbital ownership to center-fragment ownership.  Because the orbital owner has a complete periodic grid for each assigned orbital, it sends the destination core plus periodic buffer directly.  The destination retains only that core-buffer box.  After checksums and norm conservation pass, the temporary full-system orbital is released.

## 5. Localization

### 5.1 Coefficient-space objective

Compute the localization transform before full real-space materialization whenever possible.  Distributed LCFO basis blocks and coefficient batches assemble bounded 384-by-384 periodic position matrices and atomic-projector overlap matrices using the same reconstructed LCFO states and unique-core quadrature.  Contributions from every LCFO basis tail covering a core point are included before the matrix element is integrated.

Atomic s+p projections initialize and label the gauge but need not span the complete retained space; the complete periodic spread is the primary localization objective.  Within each equal-occupation block separately:

1. assemble the measured real-space LCFO metric and transform the block to a metric-orthonormal frame;
2. construct a real initial gauge from atomic projection overlap;
3. minimize the complete periodic spread;
4. use real antisymmetric generators and orthogonal exponentials;
5. map the orthogonal gauge back with exact metric bookkeeping;
6. require monotone accepted spread steps and a projected-gradient convergence gate.

The objective, gradient, and line search use the same complete periodic spread.  Fragment buffers may supply derivatives but never own integration weight.

### 5.2 Degeneracies and deterministic gauge

Degenerate or nearly degenerate LCFO eigenvectors may arrive with arbitrary real rotations.  Initialization and final ordering use projection labels, periodic centers, spreads, and a deterministic lexicographic tie-break.  A sign convention is chosen from the largest-magnitude real-space sample after localization, with ties resolved by global physical grid identifier.  Physical acceptance depends on projectors and subspaces, not individual pre-localization eigenvector identity.

Equal-occupation blocks may become small or even singleton for a metal or finite-temperature state.  Reduced localizability in that case is the physically correct consequence of preserving the density matrix: the route reports the achieved spread and rejects a failed localization gate rather than mixing unequal occupations.

## 6. Symmetry

### 6.1 What is required

The full instantaneous atomic configuration defines valid affine operations `r -> R r + t`.  Symmetry acceptance is based on:

- covariance of the LCFO Hamiltonian and occupied projector;
- invariance of the reconstructed density;
- closure of localized Wannier-center affine orbits;
- covariance of final metric, Hamiltonian, position, velocity, and nonlocal matrices.

Individual Wannier functions need not each be invariant.  They may be permuted and mixed inside a symmetry orbit.

For every accepted generator `g`, measure both `||H U_g - U_g H||` and `||P U_g - U_g P||` in the measured LCFO metric, with the boundary/interior decomposition below.  Thus a symmetric atomic configuration is not taken as proof that a finite fragment LCFO approximation retained the symmetry.  The localized basis is accepted only when it represents the same symmetry-compatible LCFO projector within the calibrated discretization error.

### 6.2 Translation subgroup and point co-group

Do not build a dense representation for every supercell translation times every point operation.  Extract the pure translation subgroup and the point co-group.  Keep fractional translations attached to point co-group representatives so screw and glide operations remain distinguishable.

Pure translations act as grid and Wannier-center permutations.  Co-group products close modulo a member of the translation subgroup; the implementation records this translation cocycle.  At Gamma the associated Bloch phase is one, while the real-space grid and Wannier-center permutation remains mandatory.  Dense orbital matrices are constructed only when required for the bounded 384-dimensional co-group action.  A common fixed point is optional metadata, never a validity requirement.

### 6.3 DC-LCFO boundary error

DC-LCFO fragment stitching may reduce pointwise smoothness at fragment faces.  Measure total, boundary-layer, and interior covariance residuals separately.  Boundary acceptance is calibrated from an independently measured LCFO face smoothness defect; the same residual in the interior is rejected.

Symmetry correction may rotate inside a fixed equal-occupation block.  It may not add orbit directions, change rank, or mix states with unequal occupations.

## 7. Ground-state reconstruction

After redistribution, assemble the DG metric and Hamiltonian from unique cores and buffer operators.  Solve the generalized eigenproblem in the fixed rank-384 basis.  The resulting occupied coefficient projector must reproduce the reference LCFO occupied projector and density.

Acceptance gates include:

- metric rank exactly 384;
- occupation spectrum and occupied rank exactly preserved (rank 128 for Si64);
- electron count;
- occupied-projector difference;
- density difference, including inversion-odd density where applicable;
- occupied-energy difference;
- center ownership and buffer coverage;
- field-off stationarity.

Only after these gates pass may V3 be published.

## 8. Checkpoint and RT

V3 records:

- LCFO basis and coefficient fingerprints;
- retained and occupied ranks;
- occupation-block and block-localization fingerprints;
- affine co-group fingerprint including fractional translations;
- translation-subgroup fingerprint;
- periodic Wannier centers and center-fragment owners;
- orbital-to-fragment redistribution receipt;
- boundary/interior covariance diagnostics;
- metric, Hamiltonian, position, velocity, and nonlocal fingerprints.

These fields define a new mandatory internal schema revision within the retained V3 checkpoint family.  A legacy V3 file lacking any mandatory provenance or acceptance receipt is rejected rather than silently reused; the production route is not renamed V4.

RT reads only an accepted V3 checkpoint.  Time propagation remains generalized-eigenvalue Exp coefficient RT.  Polarization is the primary LR/HHG observable; current is secondary.  HHG spectra are computed from polarization and shown on a semi-log scale with even-harmonic peak/dip/slope classifications.

## 9. Failure handling

Reject without fallback when:

- Gamma-real conditions are not met;
- LCFO basis rank is below the requested Wannier rank;
- the metric loses rank;
- localization mixes unequal-occupation blocks;
- occupied projector, density, or occupied energy changes beyond tolerance;
- an affine map is not grid commensurate;
- center-orbit matching is incomplete;
- redistribution loses norm, points, or buffer coverage;
- field-off RT is not stationary.

The route must not fall back to obsolete DG modes, fragment-local complement generation, displaced-structure acceptance, or current-derived HHG spectra.

## 9.1 Tolerance hierarchy

No material-specific acceptance constant is fitted to Si64.  Tolerances are derived from existing numerical controls:

- metric and projector rank: `dg_dc_metric_rank_tolerance` scaled by the measured metric condition number;
- occupied projector and density: the stricter of the final DC-SCF tolerance and propagated metric roundoff;
- interior covariance: `dg_ow_symmetry_tolerance` plus propagated eigensolver/metric error;
- boundary covariance: the independently measured LCFO face smoothness baseline plus the interior allowance;
- center matching: periodic-moment uncertainty derived from moment magnitude, localization gradient tolerance, and grid spacing;
- redistribution: collective roundoff scaled by orbital norm and message count.

Every derived allowance and its ingredients are written to evidence.  A tolerance may not be enlarged after observing whether a physical case passes.

## 10. Verification matrix

### 10.1 Synthetic and focused tests

- LCFO coefficient extraction beyond normal `nstate_tot`;
- invariance to signs and rotations within degenerate LCFO eigenspaces;
- exact preservation of integer and fractionally occupied density matrices under block-orthogonal localization;
- rejection of mixing between unequal occupations;
- orbital-batch memory bound and deterministic 1/2/4/8-rank output;
- non-divisible rank/orbital ownership and MPI count-overflow rejection;
- periodic centers crossing cell faces;
- centers outside their construction rank and final center-fragment redistribution;
- norm and density preservation across the orbital-to-spatial transpose;
- inversion, noncentrosymmetric rotations, screw, glide, `C1`, and displaced lower-symmetry structures;
- boundary-localized versus interior covariance errors;
- checkpoint write/read and provenance rejection.

### 10.2 Genuine ideal-Si64 evidence

1. strict undisplaced DC-SCF;
2. LCFO+EigenExa with at least 384 coefficient columns;
3. block-localized full-system Wannier construction;
4. redistribution to eight fragments;
5. projector, density, energy, covariance, and center-orbit gates;
6. rank-384 DG GS and V3 publication;
7. field-off RT;
8. polarization-derived linear response;
9. long-pulse Exp RT and semi-log HHG;
10. even-order peak/dip/slope classification.

Quantitative agreement is not expected from the small physical system, but qualitative symmetry and spectral morphology must be internally consistent.

## 11. Migration

Remove from the production adapter:

- fragment core-owned occupied direct sums;
- fragment-local complement construction;
- symmetry-orbit rank growth;
- provisional centers assigned by rank;
- dense enumeration of all supercell affine operations.

Retain focused low-level fixtures only where they test reusable metric, symmetry, localization, or checkpoint primitives.  Obsolete production dispatch remains forbidden.
