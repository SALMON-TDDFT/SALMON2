# WPW Density-Carrying Seed and Fixed-H Relaxation Design

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## Status and decision

This document supersedes the occupied-subspace initialization part of
`2026-07-19-wpw-fixed-hamiltonian-basis-flux-relaxation-design.md`.

The first WPW checkpoint keeps the converged DC density and every non-interface
Hamiltonian contribution fixed. It does not run a WPW density SCF. The initial
subspace is built from the core-owned projected Wannier functions selected from
the converged fragment occupied subspaces, followed by a true projection into
the nonorthogonal W+P basis.

It must not be described as a projection of the global LCFO eigenvectors.
Those vectors (`coef_wf`) are currently produced only by the Flux-inclusive
`diag_eigenexa` route that this work is intended to avoid.

## Why this source is used

The converged DC state is represented by fragment occupied subspaces. Each
fragment calculation contains buffer states and therefore exposes more
occupied eigenvectors than the number of electrons owned by its core. The
source is not the direct sum of all those eigenvectors. It is the direct sum of
projected Wannier functions whose centers are owned by the fragment core. Their
piecewise core-domain density is the density frozen for the first WPW
checkpoint.

The alternatives were rejected:

1. Running `diag_eigenexa` first would reintroduce the slow Flux-inclusive
   global eigenproblem and use a state already affected by the operator being
   introduced.
2. Treating W overlaps directly as coefficients and setting P coefficients to
   zero is not a projection in a nonorthogonal W+P basis.
3. Reconstructing global orbitals from density alone is not unique and would
   introduce an uncontrolled gauge choice.

## Core-owned Wannier source contract

For the first Si64 route, construct deterministic projected Wannier functions
without invoking the global Flux eigensolver or an external iterative MLWF
optimization. Reuse the existing bond-center/SAWF projection machinery to
build one trial per occupied Si-Si bond center, project those trials into the
converged fragment occupied subspace, and apply the polar/Löwdin factor of the
trial-overlap matrix. Reject a rank-deficient or ill-conditioned projected
trial block. The stable source ID is derived from the ordered endpoint atom
IDs and canonical periodic bond image. If endpoint ordering is reversed, also
reverse the sign of the directed image vector before canonicalization. Phase is
fixed by requiring the largest
trial-overlap component to be positive real, with stable row ID as tie-break.

Select only projected Wannier functions whose canonical bond centers lie in
that fragment's half-open core box `[lower,upper)`. Work in canonical
fractional coordinates `u=x/L-floor(x/L)`, so an exactly upper-bound center is
wrapped to the next core. Snap only values within a named center tolerance to
an exactly representable core boundary; values outside that tolerance use the
ordinary half-open comparison. Face, edge, and corner ownership is the
lexicographic result of applying the same rule in all axes. No stable source ID
may be selected by two fragments or omitted by all fragments.

The expected local source count is derived from the stoichiometric bond-center
map, not from a spatial electron-density integral and not from
`count(system%rocc>0)` over the buffered fragment calculation. For the
non-spin-polarized Si64 2x2x2 case the contract is

```text
global electrons = 256
global occupied Wannier functions = 256 / 2 = 128
core electrons per fragment = 32
core-owned occupied Wannier functions per fragment = 32 / 2 = 16
collective source count = 8 * 16 = 128
```

The value 16 is a Si64 stoichiometric and symmetry oracle derived from the
number of core-owned Si-Si bond centers. It is not inferred from integrating
the electronic density over one spatial core, which need not be an integer due
to bonding tails. The implementation must verify collectively that the selected
local counts sum to the global occupied count and that the reconstructed charge
equals the DC electron count. Nonuniform core counts are allowed when the
collective count and stable center ownership are correct. Two density checks
are distinct and both are required:

1. rotating the input occupied eigenvectors and then repeating the deterministic
   projected-Wannier construction must leave the reconstructed density and
   stable center-ID set invariant within tolerance;
2. after center ownership selection, the density reconstructed from the
   collective core-owned Wannier ensemble on the global core tiling must match
   the converged DC density within an explicit representation tolerance.

The second condition is a qualification test, not a consequence of unitarity;
Wannier tails from centers outside a core can contribute inside that core. The
selected Wannier functions must therefore be communicated to every neighboring
fragment core intersected by their buffered tails. This is an explicit
Wannier-tail halo exchange, not an owner-core-only evaluation. Each halo record
carries a stable source-Wannier ID, source fragment ID, destination core ID,
periodic image vector, and the values needed for density and W/P-overlap
accumulation. Send and receive schedules must be collectively symmetric, and a
`(source ID, destination core point, periodic image)` contribution must be
accumulated exactly once.

After halo exchange, each core reconstructs density from locally centered
Wannier functions plus received neighboring tails. The global core tiling is
compared with the converged DC density before metric projection. Overlap
contributions are routed separately by stable source and basis-row IDs. Failure
of the selected ensemble to reproduce the DC core density is fatal and must not
be hidden by renormalizing occupations. Integer
source counts are required for this first insulating, integer-occupation
production route. Metallic or fractional-occupation center assignment requires
a separate design.

### Tail support and truncation contract

Wannier functions are not assumed compact. Communicate every periodic image
available in the origin fragment buffer whose squared tail norm on a
destination core exceeds `dg_wpw_wannier_tail_norm_tolerance`. Accumulate the omitted
tail norm and omitted charge bound over deliberately discarded available
images. Values outside the finite fragment buffer cannot be bounded from the
same truncated data. Buffer sufficiency is therefore qualified independently
by both the outer-buffer-shell norm and the reconstruction residual against the
converged DC density. If either exceeds tolerance, stop and request a larger
buffer; nearest-neighbor-only fallback is forbidden.

The independent reference is a direct evaluation from the origin fragment
buffer before halo packing. A reconstruction made from the packed halo records
is compared against that reference so a missing record cannot validate itself.
The preflight must report maximum halo distance, image range, sent/received
records and values, omitted tail norm/charge bounds, and peak storage.

### Tolerance profile

Introduce named, serialized tolerances rather than hidden constants:

- `dg_wpw_wannier_center_tolerance=1d-10` in canonical fractional coordinates;
- `dg_wpw_wannier_tail_norm_tolerance=1d-12` per source/image/core record;
- `dg_wpw_wannier_density_tolerance=dg_wpw_scf_residual_tolerance` for relative
  global-core L2 density error and charge error;
- `dg_wpw_wannier_gram_cutoff=dg_wpw_metric_cutoff` for projected-trial and
  collective source Gram rank decisions.

The Si64 preflight must report sensitivity at one decade above and below the
center and tail defaults without changing center ownership, source rank, or
the accepted density within `dg_wpw_wannier_density_tolerance`. Failure of that
sensitivity check is a stop condition, not permission to tune until passing.

## Mathematical contract

Let `A=(W,P)` be the distributed W+P basis, `S=A^dagger A`, and
`Phi_src=(phi_1,...,phi_nocc)` the direct sum of core-owned projected Wannier
functions. Construct
the distributed overlap right-hand side `B=A^dagger Phi_src` from both blocks

```text
b_W = W^dagger phi_a
b_P = P^dagger phi_a
b   = (b_W,b_P).
```

The projected coefficients are the rank-qualified solution of

```text
S c = b.
```

The implementation must not replace this solve with coefficient normalization.
Let `G_C=C_raw^dagger S C_raw` and let `T` be the rank-revealing transform with
`T^dagger G_C T=I`, so `Q=C_raw T`. Because projected Wannier functions from
different buffered fragments are not assumed mutually orthogonal, the source
occupation matrix must be transformed as

```text
F_Q = T^{-1} F_src T^{-dagger}.
```

The density identity `C_raw F_src C_raw^dagger = Q F_Q Q^dagger` must pass
before fixed-H relaxation. Replacing `F_Q` by a diagonal occupation vector at
this stage is forbidden. Loss of occupied rank is fatal. Once the fixed-H
eigensolver produces its eigenbasis, physical integer occupations are diagonal
again; the resulting density change from the qualified source is recorded as
a separate fixed-H diagnostic.

Solving `S C_raw=B` and preserving density across the `C_raw -> Q` transform do
not by themselves prove that the W+P basis captured the source. Reconstruct
the physical projected density

```text
rho_projected = density[A C_raw, F_src]
rho_normalized = density[A Q, F_Q]
```

on the global core tiling. Require `rho_projected` to match both `rho_source`
and the converged DC density within `dg_wpw_wannier_density_tolerance`, including
charge. Then require `rho_normalized` to match `rho_projected` within the same
tolerance. The captured-norm deficit `abs(1-captured_norm)` must also be no
larger than that tolerance. These are publication gates, not diagnostics only.

### Distributed overlap assembly contract

The right-hand side is distributed by basis-row ownership, not by source
fragment ownership. Every fragment computes partial W and P overlaps on its
owned core points for every support basis row that is nonzero there. Those
partials are routed and summed to the canonical owner of each W/P row. A
fragment root must not assume that its source orbitals overlap only basis rows
owned by the same fragment.

The assembly oracle must include a W row owned by one root whose overlap has a
required nonzero contribution from the other fragment core. It must also
include a nonzero P block. Omitting either routed contribution must make the
oracle fail before the metric solve.

Projection diagnostics are:

- relative metric-equation residual `||S C-B||/||B||`;
- occupation-weighted captured norm `Tr(F B^dagger C_raw)` relative to the
  occupation-weighted source core norm `Tr(F Phi_src^dagger Phi_src)`, where
  `F=diag(f_a)` and `C_raw` is the converged solution of `S C_raw=B` before
  any occupied S-orthonormalization;
- occupied effective rank;
- S-orthogonality defect;
- occupation-derived charge;
- source Gram defect and condition number;
- density residual after S-orthonormalization using transformed `F_Q`;
- physical source-to-W+P projected density residual and charge error;
- captured-norm deficit with an explicit acceptance threshold;
- S-metric projector change after rotating the input occupied eigenvectors and
  repeating deterministic projected-Wannier construction.

The captured norm and metric-equation residual are properties of `C_raw`.
Rank and S-orthogonality are properties of the qualified occupied block after
normalization. Diagnostics must not substitute the normalized block into the
captured-norm formula.

The source identity oracle may apply only deterministic phase changes and
permutations to individual localized Wannier functions, transforming IDs and
`F` consistently. An arbitrary occupied unitary destroys localization and is
not required to preserve individual centers. A separate gauge oracle may rotate
the occupied fragment eigenvectors before repeating the complete projected-
Wannier construction; it compares the final center-ID set, density matrix, and
S-metric projector, not per-column phases.

### Metric-solver qualification

The distributed metric solve remains rank-local plus bounded owner/halo
exchange. It must report per-RHS residuals as well as the block residual,
detect nonfinite values and stagnation collectively, and fail at a fixed
iteration cap. The first production preconditioner is positive diagonal
Jacobi built from the owned diagonal of `S`; reject nonfinite or nonpositive
diagonal entries collectively. Report the global diagonal spread and the
per-RHS residual history as conditioning diagnostics. Introduce a block-local
preconditioner only if a focused oracle proves diagonal Jacobi insufficient,
and require a separate design review before doing so. The solver must not fall
back to a global dense solve or all-state gather.

The Si64 diagnostic now satisfies that escalation condition. With positive
diagonal Jacobi, recursive and explicitly recomputed residuals agree to about
`1e-13`, `P^dagger S P` is Hermitian to about `1e-14`, and 256 PCG iterations
do not reduce the worst RHS below `1.6e-2`. The reconstructed preconditioned
Lanczos tridiagonal has Ritz extrema `8.30e-5` and `8.81`, hence condition
number about `1.06e5`. Its relative minimum is `9.42e-6`, above the `1e-8`
metric cutoff; no cutoff-level near-null mode or RHS near-null weight was
detected. This is finite conditioning, not rank loss, an inconsistent RHS,
loss of Hermiticity, or a broken residual recurrence.

### First-stage bounded block preconditioner

The approved first-stage remedy is fragment-local coupled W+P block Jacobi,
not a claim that DG+TDDFT initial-state construction is complete. Each
canonical fragment owner forms its owned principal overlap block:

```text
S_f = [ S_WW(f)  S_WP(f) ]
      [ S_PW(f)  S_PP(f) ].
```

The block contains exactly the W and P rows canonically owned by that fragment.
Construct its principal submatrix by locating every owned W and P ID in the
corresponding required-ID array. In the current bounded layout, extract WW from
the required-W by required-W store, WP from the required-W by owned-P store,
and PP from the owned-P by required-P store. Verify the extracted lower-right
and off-diagonal blocks against Hermitian conjugacy instead of assuming that
the differently shaped stores are already an owned-by-owned matrix.

The metric Krylov vectors already contain canonical-owned rows. Consequently,
applying `S_f^{-1}` is rank-local in the current one-canonical-fragment-owner
per production rank layout; it must not add an owner exchange. Existing owner
schedules and owned/required IDs define and validate the extraction, while the
collective is used only to synchronize build/apply validation and failure
before the next distributed `S` application. A global metric gather is
forbidden. Memory and factorization cost remain bounded by one fragment block
per owner. A future layout in which one rank owns several fragments must store
and apply one independently qualified block per fragment rather than merging
them silently.

Report each block's Hermitian defect, dimension, eigenvalue extrema, and
condition number. Use a Hermitian eigendecomposition for qualification and
inverse application. An eigenvalue at or below the named local relative cutoff
is a local-rank-truncation event: fail collectively, restore the previous valid
preconditioner, and prohibit checkpoint publication. Silent truncation and
unqualified Cholesky are forbidden.

Compare diagonal Jacobi and coupled block Jacobi on the same focused MPI
fixture and Si64 preflight. Require convergence of every metric RHS below
`dg_wpw_metric_cutoff` within the fixed cap, agreement with the explicitly
recomputed residual, and at least a one-decade reduction from the diagonal
Jacobi condition estimate (`1.06e5` baseline, hence target at most `1.06e4`).
Passing this solve does not waive captured-norm, density-reproduction,
occupied-projector, or fixed-H gates.

If improvement is insufficient, stop. The next candidate is overlap-1 additive
Schwarz with neighboring-fragment support, requiring a separate ownership,
partition-of-unity, communication, rollback, and memory review. W/P-separated
block Jacobi is only a comparison oracle because it omits the coupling most
likely to control conditioning.

### Si64 block-Jacobi result

The dedicated `20260720_block_jacobi_preflight` used the corrected full-cell
buffer-shell predicate and reached the metric solve with
`outer_shell_norm=0`. Every local coupled W+P block had dimension 287,
Hermitian defect below `5.6e-17`, eigenvalue range approximately
`[1.038e-6,1.918]`, and local condition number `1.848e6`. The block-preconditioned
PCG/Lanczos estimate improved the global condition number from approximately
`1.06e5` to `2.31e3`, a factor of about 46.

This improvement is material but does not pass the strict solver target.
At the fixed 256-iteration cap the aggregate residual was `1.41e-7` and the
worst RHS residual was `3.40e-7`, above `dg_wpw_metric_cutoff=1e-10`.
Recursive and explicitly recomputed residuals agreed, the final recurrence
defect was `4.75e-13`, and no cutoff-level near-null mode was detected.

This result is not by itself a reason to replace block Jacobi. When the solve
reaches the fixed cap with finite values, positive curvature, consistent
recursive and explicit residuals, and at least one decade of improvement for
every RHS, retain the best iterate of each RHS as a **diagnostic-only metric
solution**. Continue through captured norm, projected density and charge,
S-normalization, occupied projector, zero-interface fixed-H relaxation, and
interface continuation. Mark the entire route nonpublishable: no checkpoint
or manifest may be written and no field-on or production RT may start.

Reconsider the metric acceptance rule only after those physical diagnostics
are available. If density, charge, captured norm, projector, or fixed-H
diagnostics fail, begin the overlap-1 additive Schwarz design review. If they
pass, review whether the strict `1e-10` per-RHS target is unnecessarily stronger
than the physical error contract. Do not increase the cap or relax the cutoff
silently.

The `20260720_metric_physical_diagnostic2` run reached this reconsideration
point. Its best-iterate aggregate/worst-RHS residuals were `1.39e-7` and
`2.73e-7`, but captured norm was `8.49e5`, projected-density relative error was
`1.36e6`, and projected charge was `2.717e7` instead of the fragment source
charge 32. The normalization-density residual was only `1.42e-14`, proving
that the later Löwdin/occupation transform preserved the already-broken raw
projected density. The selected source-to-DC density residual was `4.38e-2`.

Therefore the strict metric gate was not merely overconservative for this
iterate. Stop before fixed-H. Before implementing Schwarz, distinguish two
possibilities with focused diagnostics: amplification by unresolved extreme
eigenmodes of the original metric, versus an inconsistency among `S`, routed
`B`, and the real-space `A C_raw` reconstruction. Schwarz improves convergence
but cannot repair an inconsistent operator/RHS/reconstruction contract.

The `20260720_support_w_reconstruction` run resolves the latter identity for
`B`: after reconstructing `A C_raw` from all support W coefficients, routed
and direct total/W/P captures agree to `1.10e-15`. The common captured norm is
still `8.4905e5`, and the projected/normalized real-space charge is
`2.6624e14` despite S-orthogonality `5.25e-12`. This localizes the remaining
failure to the contract between the assembled metric `S` and the real-space
Gram `A^dagger A`. Split that comparison into W-W, W-P, and P-P contributions
before revisiting the iterative solver or designing overlap-1 Schwarz.

The explicit Task 2E comparison gives assembled S `8.4905e5` and real-space
Gram `8.3200e12` for the same raw coefficients. Its real-space split is WW
`8.3277e12`, WP `-7.7479e9`, and PP `4.6702e7`. Since the current WW adapter
publishes `orthonormal_ww`, the next diagnostic must test that convention
against the actual fragment-tail W fields and separate an invalid metric
assumption from W-field scaling or duplicate halo contributions.

Task 2G shows the scaling is entirely in the tail: the worst W on every
fragment has owner-core norm `1.0000` and halo-tail norm `2.5096e6`, always at
the last local W row. Canonical routing preserves total = owner-core + tail,
so this is not duplicate quadrature. The design must diagnose the buffer basis
construction/continuation and periodic packing coordinates before changing S;
an unvalidated numerical tail clamp is forbidden.

Task 2H further shows exact pre-pack/packed agreement and identifies
`transformed_spsi` as the active source. The invalid construction is therefore
upstream: core-only orthonormalized Flux-LCFO directions are continued into the
buffer without a tail-norm contract. The production correction is to honor the
fragment occupied-W design already specified above (16 projected Wannier W
rows per Si64 fragment), not to expose all 165 LCFO directions as W or clamp
their tails. P remains the augmentation space. The WW metric must be assembled
from the resulting communicated W fields unless global orthonormality is
explicitly demonstrated.

## Extra-state contract

Only occupied columns come from the DC density-carrying ensemble. Extra states
are deterministic. Before normalization, every extra column is projected out
of the occupied space:

```text
Q_extra <- Q_extra - Q_occ (Q_occ^dagger S Q_extra).
```

The extra block is then rank-qualified and S-normalized. Random replacement of
an occupied column is forbidden.

## Fixed operator and continuation

The fixed generalized eigenproblem is

```text
(H0_DC + lambda_interface H_interface) C = S C E.
```

Frozen H0 contains:

- WW kinetic, potential, and nonlocal blocks;
- WP volume and nonlocal blocks;
- PP volume and nonlocal blocks.

The complete interface block contains only:

- `ww_face_self`;
- `ww_cross_value`;
- `wp_h_face`.

Support IDs, ownership maps, halo records, and routing schedules are transport
metadata and are never scaled.

The first solve uses `lambda_interface=0`. Accepted solutions are continued to
one with rollback and step reduction. Density reconstruction is diagnostic
only and must never call `wpw_potential_step` in fixed-H mode.

## Frozen-state invariant

The snapshot must cover values, not only stored fingerprints:

- total density, Hartree, XC, and local ionic potential arrays;
- all WW component arrays and IDs;
- WP/PP volume, nonlocal, and face arrays;
- bounded H0 and interface caches;
- layout/ownership IDs and operator provenance.

Validation must be collective and shape-safe. No rank may return before a
collective entered by its peers. Allocation failure must invalidate the whole
candidate and leave no publishable checkpoint.

## Failure and publication policy

Stop before checkpoint publication if any of the following occurs:

- source occupation count or charge is inconsistent with DC provenance;
- projected-Wannier construction is not invariant to the input occupied gauge
  within tolerance;
- the selected core-owned Wannier ensemble does not reproduce the converged DC
  density on the global core tiling within the representation tolerance;
- a required neighboring Wannier tail is missing, duplicated, assigned the
  wrong periodic image, or communicated with an asymmetric halo schedule;
- omitted tail norm/charge exceeds tolerance, or the fragment buffer cannot
  establish a bounded tail error;
- a Wannier center is multiply owned, unowned, nonfinite, or ambiguously
  classified at a periodic core boundary;
- W or P overlap construction is incomplete or nonfinite;
- any support-row overlap contribution is not routed to its canonical owner;
- the metric solve has a structural failure or loses occupied rank; reaching
  the fixed cap with a qualified diagnostic-only best iterate is not a
  structural failure, but permanently prohibits publication for that run;
- any RHS stagnates, exceeds the iteration cap, or has a nonfinite residual;
- projection residual or S-orthogonality exceeds the active tolerance profile;
- source-to-W+P projected density/charge or captured-norm deficit exceeds
  `dg_wpw_wannier_density_tolerance`;
- the transformed occupation matrix is non-Hermitian/nonfinite or the density
  identity across S-orthonormalization fails;
- any frozen value, cache, ID, or provenance record changes;
- a fixed-H solve or continuation step is nonfinite or unconverged;
- `lambda_interface` has not reached one.

Long-time HHG and spectral interpretation remain out of scope.
