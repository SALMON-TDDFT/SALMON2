# Wannier90 Symmetry-Projected Line Search Design

## Problem

Wannier90 computes `doda0`, the directional derivative used by its parabolic
line search, inside `internal_search_direction`.  The SALMON symmetry patch
then projects `cdq` onto the site-symmetry-preserving tangent space.  The trial
step therefore follows the projected direction while `doda0` still describes
the unprojected direction.  The line-search model and the actual update are
inconsistent.

## Design

Keep Wannier90's optimizer and symmetry projector.  Immediately after the
site-symmetry projection, copy the projected direction back to the distributed
slice and recompute the global directional derivative from the raw gradient
and the direction that will actually be applied.  Also compute a projected
direction norm for diagnostics.

If a projected CG direction is not a descent direction, discard only its CG
history, rebuild the steepest-descent direction from the raw gradient, project
that direction, and recompute the derivative.  Do not add iteration-number
switches, material-specific step schedules, Pulay mixing, or a second optimizer.

## Verification

First add a source-route test that fails unless recomputation occurs after the
symmetry projection and before the line-search branch.  Then rebuild bundled
Wannier90 and run a short continuation from the same Si64 checkpoint with
verbose line-search output.  Confirm that the reported slope is the projected
directional derivative and compare spread monotonicity against the existing
unmodified checkpoint trajectory.  A long production run is not started until
the short diagnostic confirms the hypothesis.

