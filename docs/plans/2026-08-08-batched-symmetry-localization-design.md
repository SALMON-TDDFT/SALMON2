# Batched Symmetry-Constrained Wannier Localization Design

## Problem

The retained 48-orbital fragment block produces about 244 sparse support edges.  The current
Jacobi implementation performs a dense symmetry projection, block exponential, full-grid spread
evaluation, and collective reduction for every edge and every line-search trial.  Production stack
sampling attributes most localization time to those repeated spread evaluations.  A missing
`config.h` include also caused the production module to omit MPI collectives and run a different
localization independently on every fragment rank.

## Design

Keep the bounded symmetric strongest-neighbor graph as a localization/support diagnostic, but do
not use a non-invariant sparse graph as the variational gradient domain.  At the beginning of each
sweep, evaluate the two real localization-gradient components for all pairs in the 48-orbital local
fragment block and assemble the full anti-Hermitian negative gradient.  Apply the finite-group
Reynolds projection once to that matrix.  Orthogonal projection of the full gradient guarantees a
descent direction in the exact symmetry commutant.  Projecting only graph edges is invalid when a
dense symmetry representation maps those edges outside the selected graph; genuine Si64 exposed
this as an Armijo failure after otherwise monotone descent.

Normalize the accumulated descent generator to a bounded initial step.  Exponentiate it once per
line-search trial, apply the resulting full retained-block unitary transactionally to values,
gradients, and the published transform, and accept only a finite Armijo decrease of the collective
periodic spread.  On rejection, restore the complete backup and halve the step.  A sweep with a
gradient below tolerance converges without applying a transform.  A nonzero gradient for which no
line-search step is accepted is a hard failure rather than false convergence.

The quadratic pair work and matrices remain local to 48 orbitals rather than the 384-orbital full
system.  All fragment ranks participate in the same collective gradient and spread objective and therefore
apply the identical transform.  This preserves full-system operations that permute fragments while
allowing the Wannier data and objective contributions to remain fragment-local.

## Verification

Expose an optional spread-evaluation counter for focused tests.  A localization call may perform
the initial evaluation plus at most one line-search budget per sweep, independent of graph edge
count.  Test RED must fail against the edge-wise implementation.  Existing tests continue to gate
strict spread decrease, convergence, unitarity, exact commutation, rejection of nonconvergence, and
MPI rank independence.  Production verification repeats genuine Si64 and requires a single graph
diagnostic, monotone spread reduction, convergence, exact symmetry closure, and acceptable runtime.
