# WPW Physical-Periodic P Design

## Problem

The DC fragment orbital lives in the `domain+2B` fragment cell, while
`dc%jxyz_tot` folds that storage onto the smaller physical periodic cell. For
B=10, raw storage points separated by the 24-point physical period are not
equal because the fragment calculation itself is 32-point periodic.

## Representative rule

For every physical canonical grid point, enumerate all of its preimages in the
unwrapped fragment P. Select exactly one by the following deterministic order:

1. minimum squared grid distance to the closed fragment core box
   `[1:nxyz_domain]`;
2. lexicographically smallest unwrapped local coordinate.

The selected value is the physical representative. Populate every P preimage
of that canonical point with this value. Never average or sum images.

The helper receives the physical-cell `total_shape` and the zero-based
`fragment_origin`. For a one-based P array index `i=1:domain+2*buffer`, the
logical unwrapped grid coordinate is `g=i-buffer` (range
`1-buffer:domain+buffer`). Its one-based physical canonical coordinate is
`modulo(fragment_origin + g - 1,total_shape) + 1`, equivalently
`modulo(fragment_origin + i - buffer - 1,total_shape) + 1`, component by
component.
Only canonical points that actually have a preimage in P are grouped.  Missing
physical-cell points are valid for a partial P and do not cause failure.

For the Si64 core of 12 points and physical period 24, B=6 provides one full
physical-cell image. Increasing to B=10 adds alternative preimages but must not
change a representative already available closer to the core.

## Placement in the data flow

Apply physical periodization to the assembled, reordered occupied-orbital P
array immediately after `unwrapped_buffer_order`. Perform projected-W
transformation and analytic finite-difference gradients only after this step.
Because both are linear operations, values and gradients then share the same
physical-periodic P contract.

The raw DC fragment orbitals remain unchanged; their fragment-cell boundary
condition is not overwritten. Only the WPW bootstrap copy is periodized.

## Invariants

- Core values are unchanged: a core point is always distance zero and wins.
- Supported geometry requires `nxyz_domain <= total_shape` componentwise, so
  distinct core points cannot be physical-period aliases; reject violations.
- B=6 and B=10 produce identical canonical representatives wherever their raw
  nearest-core values agree.
- Every physical-period alias in P is exactly equal after construction.
- Periodization is idempotent.
- No values are summed. Canonical points absent from a partial P are allowed.
- D remains P eroded by the derivative radius; no new buffer is requested.

## Verification and stopping

Unit tests cover B=5/B=6/B=10, absolute origin/canonical seam mappings,
multi-axis aliases, lexicographic tie breaking, core preservation, idempotence, unsupported
geometry, partial coverage, and non-finite input. Production adds a
`physical_periodic_p` collective stage. For every B>4 it checks equality of
all aliases that are present in P. Full canonical-cell extraction and the
spread diagnostic are a separate postcondition and run only when P covers the
whole physical cell. Fresh B=6/B=10 runs must produce the same 128-W spread
distribution within the existing named tolerances. Existing outer-shell,
fixed-H, publication, and RT gates remain unchanged.
