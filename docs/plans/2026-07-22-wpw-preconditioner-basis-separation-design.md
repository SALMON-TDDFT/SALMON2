# WPW Preconditioner versus Basis-conditioning Design

## Goal

Separate failure caused by the diagonal eigensolver procedure from failure
caused by an intrinsically ill-conditioned W/P coordinate basis.

## Design

Keep the current preconditioned production route as the default.  Add an
explicit diagnostic input switch, defaulting to enabled, that omits the
optional `wpw_precondition` callback only when the user requests a comparison
run.  The unpreconditioned route uses the identical normalized occupied-W
descriptor, W/P metric and Hamiltonian, deterministic extra states, tolerance,
iteration cap, and fixed-H gates.  It remains nonpublishable unless every
existing production qualification passes.

For the enabled route, instrument the diagonal preconditioner at the same
selected window iterations used by the state-residual diagnostic.  For each
raw-worst occupied and extra state, report separately for W and P rows:

- minimum and maximum `abs(H_ii - epsilon S_ii)`;
- number of positive, negative, and floor-clamped denominators;
- the scale-normalized minimum denominator;
- maximum absolute inverse multiplier actually applied.

The preconditioner callback receives no inner-iteration number, so the
selected diagnostic state IDs will be passed through a small diagnostic state
owned by the matrix-free driver immediately before the callback and cleared
afterward.  The production callback may read it only for logging; numerical
output remains bitwise identical when diagnostics are enabled.

## Interpretation

- If disabling the preconditioner materially improves both occupied and extra
  residual histories, the procedure is the immediate defect.
- If both routes stall similarly while the metric is strongly ill-conditioned,
  basis conditioning or the LOBPCG search-space update is primary.
- If only extra states improve, revise extra-state completion before changing
  occupied acceptance.
- Floor hits or very small normalized denominators correlated with large
  amplification justify a shifted/bounded preconditioner in a later design.

No tolerance relaxation, iteration increase, basis rotation, or checkpoint
policy change is part of this diagnostic task.

## Verification

Use TDD for input default/broadcast/validation, callback omission, diagnostic
state transaction, and denominator classification.  Run focused MPI tests,
full build, code review, then matched preconditioned and unpreconditioned B=6
runs in fresh directories.
