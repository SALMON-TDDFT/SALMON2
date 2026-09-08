# WPW Tail Diagnostic Design

## Decision

The prepared-P outer-shell norm is a convergence diagnostic, not a mandatory
physical gate. Projected Wannier columns are already normalized on unique
canonical points. Continue descriptor construction when the outer ratio is
structurally valid under the rules below, even if it exceeds the configured
reference tolerance.

Keep the current total norm, outer-shell norm, ratio, and tolerance in the
log, and add an explicit `status=warning` when the ratio exceeds tolerance.

Structural validity is separate from that warning:

- total norm must be finite and strictly positive;
- outer-shell norm must be finite and non-negative;
- outer-shell norm must not exceed total norm, allowing
  `100*epsilon*max(1,total)` for summation roundoff;
- the reported ratio is clamped to `[0,1]` only after those checks.

Thus zero outer norm is valid and ideal. Non-finite or structurally impossible
values and MPI failure remain collectively fatal. Normalization failure,
metric/Gram failure, density/electron-count failure, and all later physical
gates remain fatal. No tolerance is changed.

This is intentionally approximate: B controls the prepared support, while
tail diagnostics provide evidence for later convergence studies. If
normalized-W density, metric, electron count, or fixed-H calculations fail,
tail treatment and localization are reconsidered then.
