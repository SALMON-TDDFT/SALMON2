# WPW Physical-Cell Wannier Normalization Design

## Problem

For B=10, raw occupied fragment orbitals are normalized in a 32-point
fragment cell. Physical periodization selects one representative for each
point of the 24-point physical cell. The resulting projected Wannier norms
are 0.987447--0.990242, so the unchanged production requirement
`abs(N_w-1) <= 1e-8` correctly rejects the descriptor.

## Decision

Normalize every projected Wannier independently on the physical cell. After
the first projected-P construction and canonical extraction, compute
`scale_w=1/sqrt(N_w)`. Reject non-positive or non-finite norms. Do not mix W
columns.

Apply the same scale exactly once to:

- the corresponding column of `polar_transform`;
- the already constructed fragment core values;
- the already constructed projected P values.

Analytic gradients are constructed later from the scaled `polar_transform`,
so values and gradients receive the identical factor. Rebuild canonical
links after scaling and retain the existing `abs(N_w-1) <= 1e-8` check.

## Invariants

- Stable W IDs, centers, and column ordering do not change.
- Raw `spsi` and physical-period representative selection do not change.
- No unitary/localization rotation is introduced.
- Core values, P values, and analytic gradients use one column scale.
- B remains a user convergence parameter and no larger buffer is requested.
- The outer-shell gate and its tolerance remain unchanged.

## Verification

Unit tests verify exact complex column scaling, invalid norm rejection, and
idempotent unit-norm input. Source-contract tests require scaling after the
first physical-cell norm measurement and before analytic gradients and
descriptor publication. Fresh B=6 and B=10 runs must pass canonical extraction
and the unit-norm spread gate; their normalized spreads are then compared as a
buffer-convergence result rather than assumed identical a priori.
