# Point-Cogroup Covariant Center Gauge Design

## Goal

Construct a deterministic post-Wannier90 gauge whose translation-character
sectors, Gamma pairing, point-cogroup cocycle, and periodic Wannier centers are
simultaneously covariant, without introducing replicated full-rank matrices.

## Root cause

The Si64 post-character basis realizes affine operation 2 with unitarity defect
`1.39e-14`, but its elementwise monomial defect is `8.76e-1` and all weight
leaks outside tolerance-compatible center blocks.  The current gauge therefore
preserves the retained subspace while mixing orbitals belonging to different
Wannier centers.  Tightening or relaxing the center matcher cannot repair this
missing point-center gauge.

## Mathematical construction

Let each translation character sector have multiplicity `m`.  In every sector
construct the three projected periodic-position matrices

```text
Z_a(c) = C_c^H exp(i 2 pi r_a) C_c,  a=1,2,3.
```

The point representative `p` maps character `c` to `p c` and acts internally
through the already measured unitary intertwiner `D_p(c)`.  Covariantize the
position tuple by streaming the 48 point representatives and averaging each
transported component with its integer Cartesian rotation.  The resulting
three-matrix tuple transforms exactly under the supplied point-cogroup action.

Choose a deterministic internal frame from this tuple.  Diagonalize a bounded
Hermitian discriminator formed from fixed irrational coefficients multiplying
the Hermitian real and imaginary parts of the three matrices.  Detect every
tolerance-degenerate block.  Within such a block, project and diagonalize the
existing LCFO physical operator.  If that remains degenerate, retain the whole
block as an unresolved physical multiplet and transport one common gauge to
all related characters; do not accept independently chosen LAPACK gauges.

For a non-self-conjugate character pair, construct one side and generate its
partner through the measured Gamma sewing.  For a self-conjugate sector, reuse
the existing antiunitary fixed-frame construction.  The inverse translation
transform is then applied exactly once per character.

## Architecture and memory

Add one primitive that operates on one character orbit at a time.  Inputs are
row-owned sector coefficients, the three local periodic phase arrays, compact
point-cogroup character/action metadata, the LCFO projected operator, Gamma
receipts, and catalog provenance.  Outputs are the rotated row-owned sector,
its conjugate partner when required, and diagnostic receipts.

Persistent large data remain the existing row-owned sector and spatial basis.
New dense storage is `O(m^2)` per periodic component and point-operation tile;
for Si64, `m=12`.  Point representations and projected position matrices are
streamed.  No replicated `Nwann x Nwann`, all-sector tensor, or `Ngrid x T`
array is added.

## Determinism and degeneracy

Canonical character ordering and Task1 catalog fingerprints define traversal
order.  Eigenvalue clusters use adjacent-gap chaining and rank-agreed
boundaries.  Nondegenerate eigenvectors receive the existing global-row phase
fix.  Exact multiplets are never split by numerical noise: their common
projector is fingerprinted and the same transported block gauge is used across
the entire point-character orbit.

## Validation and receipts

Before accepting the new gauge, measure collectively:

- sector Gram-I and internal rotation unitarity;
- periodic-position covariance under every point generator;
- Gamma fixed/pair residual;
- factored point-cogroup cocycle closure;
- inverse-transform imaginary defect;
- affine center-orbit residual;
- pre/post total spread and maximum spread increase;
- maximum internal rotation norm;
- row-owned and `O(m^2)` workspace peaks;
- a provenance fingerprint binding catalog, position tuple, LCFO operator,
  Gamma sewing, transported multiplet projectors, and final center projector.

The existing center tolerance is unchanged.  Nonfinite input, rank disagreement,
unsafe extent arithmetic, allocation failure, unresolved nonphysical splitting,
or failure to improve center closure rejects collectively.

## Tests

Focused MPI fixtures cover:

- Z4 with a non-self-inverse character and directed point action;
- Z2 x Z2 with self-conjugate characters;
- an internal unitary mixing two distinct centers, repaired to zero leakage;
- an arbitrary rotation inside an exactly repeated-center multiplet;
- LCFO secondary discrimination of a periodic-position degeneracy;
- a deliberately unresolved/inconsistent block that rejects;
- reversed point/character action orientation;
- Gamma-pair and self-conjugate reality;
- duplicate/missing ownership, rank-disagreeing metadata, nonfinite and
  finite-huge payloads, and checked workspace boundaries;
- identical physical receipts on MPI 1/2/4/8.

Production integration processes one point-character orbit at a time before
the inverse translation accumulator finalizes.  The existing center-gauge
diagnostic remains as an independent postcondition and must report leakage at
or below tolerance before publication.
