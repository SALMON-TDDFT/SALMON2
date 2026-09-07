# Wannier90 Symmetry-Constrained Monotonic Backtracking Design

## Motivation

For the Si64 Gamma-only symmetry-constrained route, Wannier90's analytic
directional derivative ceases to track the finite-difference spread derivative
near the localized basin.  The mismatch is state dependent, so a fixed
normalization correction is not valid.  The existing parabolic line search can
therefore accept a step that increases the actual spread.

## Algorithm

Use derivative-free monotonic backtracking only when site symmetry is active
and a fixed step was not explicitly requested.  Preserve the current `U` and
`M` through Wannier90's existing `u0_loc` and `m0_loc`/scratch-file storage.
Start from the last accepted trial step, evaluate the real spread, and accept
when it does not exceed the starting spread beyond roundoff.  Otherwise restore
the saved state, halve the step, and retry.  Stop reducing at a scale derived
from machine precision rather than an iteration-count or material-specific
threshold.  If no decreasing step exists, restore the original state, take a
zero step, and reset CG history.

The accepted step becomes the next iteration's initial step.  It is never
increased automatically: the initial large step cools quickly, while rejected
steps produce gradual, data-driven attenuation near the constrained minimum.

## Scope and memory

The unconstrained upstream Wannier90 line search and explicit `fixed_step`
behavior remain unchanged.  No new orbital-sized array is allocated.  Temporary
diagnostics used to investigate slope normalization and anti-Hermiticity are
removed from the production path.

## Verification

Source-route tests require symmetry-only backtracking, rollback before each
retry, monotonic acceptance, machine-precision termination, CG reset on a zero
step, and absence of the temporary O(N²) diagnostic loop.  A short Si64 restart
must show non-increasing spread at every accepted iteration before any long run
is started.

