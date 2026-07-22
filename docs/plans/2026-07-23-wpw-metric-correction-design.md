# WPW Metric-aware Correction Design

## Purpose

Test whether the fixed-H plateau is caused by using residuals as correction
vectors without accounting for the non-identity W/P metric.

## Chosen comparison

Reuse the already initialized fragment-local metric block eigendecomposition and
form the LOBPCG correction

```text
z = S_block^{-1} r.
```

`S_block` is the owned W/P principal block used by the occupied-seed metric
solver.  Its initialization validates Hermiticity, positive definiteness,
cutoff, operator epoch, layout fingerprint, dimensions, and finite spectral
data.  Applying its inverse requires no additional H/S operator call and no
nested iteration.

Add a default-off `yn_dg_wpw_metric_preconditioner` control.  It is mutually
exclusive with the existing diagonal `yn_dg_wpw_preconditioner`; enabling both
is an input error.  The default input therefore preserves current behavior.
The first comparison explicitly sets diagonal preconditioning off and metric
block preconditioning on.

## Alternatives

- A global iterative `S^{-1}r` would be more accurate but introduces a nested
  160-RHS metric solve at every inner iteration and inherits its strict
  convergence problem.
- A state-dependent shifted `(H-lambda S+sigma S)^{-1}r` correction is closer to
  the ideal correction equation, but requires a justified shift and stable local
  H blocks.  It is the next option only if the SPD metric correction is
  insufficient.

## Invariants and interpretation

Search history, basis, cutoff, tolerance, 160-inner cap, convergence predicate,
fixed-H/continuation transaction, and publication gates remain unchanged.  The
existing selected state residual diagnostics measure the correction norm ratio;
search metric-mode and Ritz-consistency diagnostics remain active.

Improvement over Task 16 identifies the missing S-metric correction geometry as
an important cause.  No improvement or deterioration rules out pure metric
scaling and motivates the shifted state-dependent correction equation.  Neither
outcome by itself authorizes relaxing residual tolerances or publication gates.

## Verification

Use TDD for input plumbing, mutual exclusion, callback selection, and runtime
application.  Extend the MPI fixture with a non-identity SPD metric whose exact
inverse correction is known.  Run focused tests, the full build, code review,
then a fresh B=6 comparison with metric block correction enabled.
