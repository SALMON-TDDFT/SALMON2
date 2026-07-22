# WPW Window State-residual Diagnostic Design

## Context

The normalized B=6 occupied-W route now reaches the zero-interface fixed-H
window solve.  Metric orthogonality remains near `1E-14`, but the maximum
residual over all 160 retained states stalls near `1.7E-3` and reaches the
160-iteration limit (`info=40`).  The current scalar maximum cannot determine
whether the physically occupied 128-state projector is stalled or whether the
32 deterministic extra states alone control the failure.

## Decision

Keep the solver, tolerance, iteration cap, acceptance rule, and checkpoint
boundary unchanged.  Extend `solve_window` with the occupied-state count and
emit a diagnostic summary at iterations 1, 2, 4, 8, every 16 iterations, and
the final iteration:

- maximum residual and worst one-based state index in the occupied block;
- maximum residual and worst one-based state index in the extra block;
- norm of the preconditioned residual for each block and its ratio to the raw
  residual norm at the block's raw-worst state.

Residual norms use the same production global Gram operation as the existing
maximum.  Preconditioned norms are diagnostics only and are computed only at
the selected iterations to limit collective cost.  Invalid/nonfinite
diagnostic values remain fatal because they expose an invalid solver state.
No convergence decision may consume these new split values in this task.

## Alternatives considered

1. Relax the all-state tolerance or accept the occupied projector directly.
   Rejected until the split residual proves which block controls failure.
2. Increase the 160-iteration cap.  Rejected because the observed plateau may
   be a seed/preconditioner defect rather than insufficient iterations.
3. Dump all 160 per-state values every iteration.  Rejected as excessive; the
   block maxima, stable state IDs, and selected-iteration preconditioner ratios
   provide the decision evidence with bounded output.

## Verification

Add a pure summary helper with tests for occupied-worst, extra-worst, ties,
invalid partition, and nonfinite norms.  Add source-contract checks that the
diagnostic cannot alter convergence.  Run matrix-free/fixed-H focused tests,
the MPI build, code review, and a fresh B=6 physical run.
