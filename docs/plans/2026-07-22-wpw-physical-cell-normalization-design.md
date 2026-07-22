# WPW Physical-Cell Wannier Normalization Design

## Problem

For B=10, raw occupied fragment orbitals are normalized in a 32-point
fragment cell. Physical periodization selects one representative for each
point of the 24-point physical cell. The resulting projected Wannier norms
are 0.987447--0.990242, so the unchanged production requirement
`abs(N_w-1) <= 1e-8` correctly rejects the descriptor.

## Decision

Normalize every projected Wannier independently using the canonical points
available in prepared P. Count each physical-period canonical point at most
once and multiply the sum by the real-space quadrature weight `system%hvol`.
For B large enough to cover the physical cell (B>=6 in Si64), this is
the complete physical-cell norm. For a partial P such as B=5, missing
canonical points are the truncated tail and contribute zero; this is the
finite-buffer approximation whose convergence is controlled by B. Never
request an additional point. After the first projected-P construction,
compute `scale_w=1/sqrt(N_w)`. Reject non-positive or non-finite norms. Do not
mix W columns.

Apply the same scale exactly once to:

- the corresponding column of `polar_transform`;
- the already constructed fragment core values;
- the already constructed projected P values.

Analytic gradients are constructed later from the scaled `polar_transform`,
so values and gradients receive the identical factor. For full-cell P,
rebuild canonical links after scaling and retain the existing
`abs(N_w-1) <= 1e-8` spread check. Partial P cannot provide the spread links,
but its occupied-W descriptor is still normalized and constructed.

The owner fragment root places its first-pass norms into the existing global
W vector. `dc%icomm_tot` collects that vector, and every rank extracts the
same `[source_offset+1:source_offset+source_count]` slice. All ranks validate
the slice collectively before any mutation and then apply identical scales.

After scaling, update the local slice of `source_values`, rebuild
`occupied_w_p`, and clear canonical values plus all norm/link local, partial,
and global arrays before the second pass. Reassemble and globally collect the
unique-canonical norm for every B, including partial B=5, and require
`abs(N_w-1)<=1e-8`. Full-cell P additionally rebuilds links and spread. This
prevents first-pass data from being retained or accumulated twice.

## Invariants

- Stable W IDs, centers, and column ordering do not change.
- Raw `spsi` and physical-period representative selection do not change.
- No unitary/localization rotation is introduced.
- Core values, P values, and analytic gradients use one column scale.
- B>4 remains a user convergence parameter and no larger buffer is requested.
- The outer-shell gate and its tolerance remain unchanged.

## Verification

Unit tests verify exact complex column scaling, invalid norm rejection,
idempotent unit-norm input, and non-unit-weight unique canonical-point norm
assembly for both partial and aliased P. Source-contract tests require total-communicator norm
collection and validation before scaling, explicit second-pass clearing, and
scaling before analytic gradients and descriptor publication. Fresh B=6 and
B=10 runs must pass canonical extraction and the unit-norm spread gate; their
normalized spreads are then compared as a buffer-convergence result rather
than assumed identical a priori.
