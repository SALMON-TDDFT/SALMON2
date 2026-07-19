# WPW Nonlocal Projector Blocks Design

## Problem

The production WPW operator imports a fragment-local WW nonlocal block, but
its WP and PP volume quadrature contains only kinetic and local-potential
terms.  It also cannot represent projector-mediated WW coupling beyond the
fragment-local block.  This violates the approved separable pseudopotential
contract and permits near-redundant W/P combinations to see inconsistent
Hamiltonians.  The Si64 Ritz trace consequently develops very deep negative
states and fails to converge despite a well-conditioned retained metric.

## Mathematical representation

For every global projector `q`, use SALMON's existing convention

```text
u(q,a) = hvol sum_grid conj(beta_q(grid)) B_a(grid)
H_nl(a,b) = sum_q conj(u(q,a)) rinv_uvu(q) u(q,b).
```

The projector grid values and `rinv_uvu` come from the existing pseudopotential
grid.  No new radial interpolation, normalization, or phase convention is
introduced.

## Ownership and locality

Each fragment root constructs sparse projector-overlap records only for basis
IDs it owns:

- W is piecewise-DG and is integrated only on its owned fragment core.
- P is evaluated analytically over its bounded window support.
- zero projector overlaps are omitted.

Owned overlap records are sent only to the existing bounded support peers.
Each row owner then intersects projector IDs and forms all required WW, WP,
and PP row/support-column blocks.  Communication degree remains bounded by
the existing support stencil; no global basis vector, dense global matrix, or
all-fragment candidate gather is allowed.

The new WW nonlocal outer product becomes authoritative.  During GS
construction, its owned-owned subblock is compared against the existing LCFO
fragment-local nonlocal matrix and fails closed on a convention mismatch.
The existing LCFO matrix is not added again.

## Lifetime

Nonlocal blocks are fixed for field-off GS and are stored separately from
kinetic, local-potential, and face contributions.  SCF potential updates
replace only the local-potential volume terms and re-add the fixed nonlocal
arrays.  Checkpoint H blocks therefore contain the complete operator without
requiring pseudopotential reconstruction during RT handoff.

## Failure behavior

Reject nonfinite overlaps, duplicate `(basis,projector)` records, inconsistent
projector weights, missing support records, LCFO local-block disagreement,
and any collective routing failure.  Diagnose the detecting rank before the
collective fail-closed path.

## Verification

1. A dense toy oracle with two W, two P, and three signed projector weights
   must reproduce every WW/WP/PW/PP entry.
2. A two-rank fixture must route sparse overlap records and reproduce the same
   dense oracle without global basis storage.
3. Potential replacement must preserve fixed nonlocal blocks exactly.
4. Existing bounded H/S and checkpoint fixtures must remain unchanged.
5. A fresh Si64 preflight must eliminate the nonlocal inconsistency and the
   unphysical Ritz collapse before production GS is permitted.
