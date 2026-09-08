# Boundary-Smooth Global LCFO Pencil Design

## Correction to the occupied-projector design

The DC-LCFO occupied orbitals are eigenstates of fragment Hamiltonians whose
density and effective potential are not smooth across fragment boundaries.
They therefore contain both the physical full-system symmetry and artificial
symmetry introduced by the fragment tiling, buffer truncation, and boundary
conditions.  Group averaging the LCFO occupied projector cannot remove an
artifact which already commutes with the physical group.

Consequently LCFO orbitals must not define the physical occupied projector.
They are used only as local expansion data in a larger overlapping
core-plus-buffer basis.  The global occupied space is obtained from a new,
boundary-smooth generalized eigenproblem.

## Construction

For every fragment, define a partition-of-unity weight on its core plus buffer.
The weights are one in the uncontested core, taper smoothly through the overlap,
vanish at the outer buffer boundary, and sum to one at every global grid point.
The taper and its first discrete derivative must match across paired fragment
faces.  Fragment identity and rank may not enter the physical value of a
weight.

Use the normal DC-LCFO coefficients and all buffer-composed complete local
`s+p` channels to evaluate local basis functions.  Accumulate, without
materializing full-system orbitals,

`S_ab = sum_f <sqrt(w_f) phi_fa | sqrt(w_f) phi_fb>`

and the correspondingly stitched Hamiltonian and density matrices.  Kinetic
terms include derivatives of the weighted functions; omitting the weight
gradient would reintroduce a boundary delta-like defect.  Local, nonlocal
pseudopotential, Hartree, exchange-correlation, and overlap contributions use
the same unique physical-grid ownership and buffer coverage rules.

The matrices are then averaged under the complete atomic space-group action:

`A_sym = |G|^-1 sum_g D_g^dagger A D_g`, for `A in {H,S,rho}`.

The operation representations are streamed by generators and checked against
the full affine product/cocycle catalog.  Symmetrization is performed on
distributed matrix tiles; no full-system wavefunction tensor and no
`|G|*Nbasis**2` tensor is retained.

Solve

`H_sym C = S_sym C epsilon`

with the existing distributed generalized algebra and EigenExa.  Select the
lowest 128 occupied states by complete degenerate/irreducible blocks.  Reject
a Fermi boundary that cuts a block, rank loss, an indefinite overlap, or a
non-real Gamma gauge.  These states, not the fragment LCFO eigenstates, form
the density-carrying part of the Wannier seed.  Direct-sum them with the
complete 256-channel `s+p` complement and require total rank 384.

## What symmetry can and cannot repair

Space-group averaging removes components inconsistent with the atomic
symmetry.  It cannot diagnose a boundary artifact that is itself symmetric.
Therefore acceptance also requires direct smoothness and reconstruction
receipts before and after matrix symmetrization:

- partition sum and discrete-gradient continuity;
- density and density-gradient jumps across every paired core face;
- kinetic-energy boundary contribution;
- electron-count and occupied-density change relative to normal DC;
- interior and buffer-shell density differences;
- commutators of `H`, `S`, and `rho` with every affine generator; and
- occupied residual `||H C - S C epsilon||`.

The boundary receipts are physical gates, not tolerances inferred from the
same discontinuous LCFO solution.  Small Si64 is used for qualitative symmetry
and harmonic-selection validation only; it is not a claim of bulk quantitative
accuracy.

## Memory and parallelism

Spatial ranks retain fragment-plus-buffer basis tiles and row-owned matrix
tiles.  Contributions are reduced to their matrix-tile owners and released
before the next orbital tile.  EigenExa consumes its padded local cyclic
blocks directly.  Peak current bytes are measured on every rank and reduced;
estimated or zero workspace receipts are invalid.

The preprocessing occurs once before RT.  Its memory bound is independent of
the number of affine operations and does not replicate the real-space basis.

## Failure and fallback policy

Reject rather than:

- using raw LCFO occupied eigenstates as the physical projector;
- excusing boundary jumps through a measured LCFO allowance;
- dropping weight-gradient kinetic terms;
- polar-unitarizing a leaky representation;
- splitting a symmetry-degenerate occupied block; or
- falling back to unconstrained Wannier90 or the removed experimental routes.

Normal DC LCFO+EigenExa remains unchanged.  Only the accepted overlapping-
Wannier preprocessing consumes the stitched global pencil.  Downstream remains
converged symmetry-adapted Wannier90, V3, and generalized-eigenvalue Exp-only
coefficient RT.

## Acceptance

Every task uses a genuine RED, focused MPI 1/2/4/8 verification, specification
and code-quality reviews, resolution of all Critical/Important findings, and a
clean-first committed-parent prerequisite overlay.  After ideal undisplaced
Si64 GS/V3 passes, resume polarization-primary LR and long-pulse HHG.  The
semilog polarization spectrum must explicitly determine whether H2 and H4 are
peaks or dips; current remains secondary.
