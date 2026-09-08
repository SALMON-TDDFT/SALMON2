# DG Fragment-PW Orthonormalization Design

## Goal

Make the mixed fragment + plane-wave basis explicitly orthonormal between the fragment (`F`) and plane-wave (`P`) subspaces at startup and after adaptive basis updates, so runtime propagation no longer needs to carry dense `FP/PF` overlap couplings as a first-class path.

## Current Problem

The current mixed-basis path computes and stores `S_mat_frag_pw = <F|P>` and propagates with nonzero `FP/PF` overlap couplings. A QR-based rank filter removes obviously redundant PW columns, but it does not orthogonalize the surviving PWs against the fragment basis. As a result:

- mixed overlap/operator application must keep `FP/PF` terms alive,
- overlap solves remain more complicated than the fragment-only case,
- basis-update paths preserve the mixed representation in-place instead of rebuilding an orthogonal mixed basis,
- the code already hints at a better direction via `pw_orthogonalized`, but that storage is currently unused.

## Recommended Approach

Use explicit fragment-projected orthonormalization of the PW subspace:

1. Build `S_FF` and `S_FP`.
2. Solve for the fragment projection coefficients `C = S_FF^{-1} S_FP`.
3. Form the fragment-orthogonalized PW basis
   `P_perp = P - F C`.
4. Build the PW self-overlap in the orthogonal complement,
   `S_perp = <P_perp|P_perp>`.
5. Orthonormalize `P_perp` internally, keeping only numerically independent directions.
6. Rebuild `S/H` couplings in this transformed basis.

With this representation, `S_FP` should become numerically negligible, `S_PP` should be close to identity, and propagation can be simplified toward block-local fragment evolution plus diagonal/identity PW handling.

## Scope

This design covers:

- startup mixed-basis preparation,
- mixed-basis diagonalization,
- adaptive basis update refresh for mixed runs,
- coefficient transfer between old and new mixed bases.

It intentionally does not yet remove all dense fallback codepaths. Those can be deleted after the orthogonalized path is validated.

## Data Model Changes

Add persistent storage for the orthogonalized PW representation in `s_dg_fragment_rt`:

- a real-space orthogonalized PW basis array replacing the current unused placeholder,
- optional transformation matrices from raw PW coefficients to orthogonalized PW coefficients,
- diagnostics for post-orthogonalization overlap quality.

The existing `S_mat_frag_pw` and `H_mat_frag_pw` can remain during the transition, but after orthogonalization they should represent the transformed basis and be numerically near-zero in overlap.

## Algorithm Details

### Startup

After `init_plane_wave_basis`, construct the raw PW basis as today, then immediately orthogonalize it against the fragment basis. All mixed operators (`S_mat_frag_pw`, `H_mat_frag_pw`, PW diagonal terms, and optional startup diagonalization) must be built in this orthogonalized basis.

### Basis Update

When the fragment basis changes, the old orthogonalized PW basis is no longer guaranteed to remain orthogonal to the new fragment basis. Re-run the same orthogonalization pipeline after the fragment basis update and before coefficient projection.

### Coefficient Transfer

Coefficient projection must use the updated orthogonalized mixed basis. In particular, the PW part of `U_new` and the mixed overlap used in `project_coefficients_mixed_state` must be consistent with the transformed PW basis, not the raw plane waves.

## Runtime Consequences

If implemented correctly:

- `S_FP` in propagation should be removable or reduced to diagnostics/fallback,
- restricted local overlap solves become much more faithful in mixed runs,
- `FP/PF` momentum and overlap special cases can be reevaluated and likely simplified,
- scalability improves because runtime coupling between fragment and PW subspaces is pushed into infrequent rebuild steps instead of every propagation step.

## Validation

Required checks:

- `max |S_FP|` after orthogonalization is near machine precision or the configured tolerance,
- `S_PP` is close to identity for the kept PW subspace,
- mixed startup and basis-update paths preserve norm and do not explode the condition number,
- RT propagation still advances through `after-overlap-solve` in mixed runs,
- energies and occupied-state norms remain stable relative to the current baseline.

## Risks

- Real-space storage of orthogonalized PWs may increase memory if done naively.
- Basis-update logic may silently mismatch old/new PW coefficient meaning unless transformation ownership is handled carefully.
- If fragment overlap `S_FF` is ill-conditioned, orthogonalization against `F` must reuse the same regularization policy already used elsewhere in the code.

## Chosen Direction

Implement true `F-P` orthonormalization at startup and basis-update time, and then migrate propagation code to assume the mixed basis is already orthogonalized. Keep the current `FP/PF` runtime path only as a temporary fallback until validation is complete.
