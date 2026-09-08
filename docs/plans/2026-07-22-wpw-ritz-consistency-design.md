# WPW Ritz Consistency Diagnostic Design

## Purpose

Determine whether the fixed-H residual plateau is introduced at the reduced
Rayleigh--Ritz update boundary, without changing the basis, search history,
preconditioner choice, cutoff, tolerance, or convergence/publication gates.

## Chosen approach

Keep the current history-retaining, no-preconditioner solver and add observation
at selected inner iterations.  Immediately after the reduced solve, form the
updated Ritz images from already available arrays,

```text
Q_new  = Z C
HQ_new = HZ C
SQ_new = SZ C
R_new  = HQ_new - SQ_new diag(lambda)
```

This adds no Hamiltonian or metric application.  Measure global per-state norms
of `R_new`, the maximum occupied/extra residual, and `Q_new^H SQ_new-I`.
At the beginning of the next inner iteration, the production code applies H and
S directly to the stored `Q_new`; form the corresponding pre-diagonalization
residual using the carried eigenvalues.  Store the complete rank-local
post-update residual vector, then use the production global Gram to measure the
per-state norm of `R_direct-R_post`.  Report occupied/extra maxima and stable
worst-state indices.  A uniform pending flag is set and cleared on every rank.
Define each relative vector defect symmetrically as
`||R_direct-R_post||/max(||R_direct||,||R_post||)`; when both norms are zero the
relative defect is zero.  Block maxima use the first state index on ties.

The selected schedule is defined by the direct comparison iteration.  To report
comparisons at inner 32, 96, and 160, arm the pending transaction after updates
31, 95, and 159 respectively.  Never compare post-update 160 with direct 160;
they represent different Ritz states, and post-update 160 has no following
iteration under the current cap.

For the same updated Ritz vectors, explicitly construct the reduced-coordinate
residual

```text
R_red = H_red C - S_red C diag(lambda)
```

and summarize it per state and per occupied/extra block.  This gives a matched
comparison between projected Ritz accuracy and the physical Ritz residual,
rather than comparing an unmatched global scalar with block maxima.

## Interpretation

- A large post-update physical residual with a small reduced-coordinate residual
  means the reduced solve is accurate only inside the projected space; this is
  expected Ritz behavior and shifts attention to correction-direction quality.
- A large same-state vector defect between post-update `HZ*C`/`SZ*C` residual
  and the next direct H/S residual identifies an operator or update recurrence
  inconsistency.  Similar block maxima alone are not treated as evidence.
- Agreement at that boundary, followed by weak reduction after the next
  Rayleigh--Ritz step, rules out stale Ritz images and points to the available
  correction subspace rather than the occupied-W representation itself.

The diagnostic is observational and runs only after the existing convergence
test has failed.  It adds no H/S/preconditioner call.  Before any diagnostic
Gram, local finite/shape validity is reduced collectively; no rank may return or
skip a collective independently.  Pending state is cleared uniformly on success
and error.  Non-finite diagnostic values remain fatal, but no diagnostic
magnitude is used as a convergence or publication gate.

## Alternatives rejected

Reapplying H and S immediately after every update would test the same boundary
but duplicate the expensive production operation.  Logging every state at every
iteration would create excessive output; selected iterations and occupied/extra
maxima are sufficient for the current hypothesis.

## Verification

Add source-contract checks for the two boundary measurements and unchanged
convergence predicate.  Extend the deterministic MPI fixture to exercise the
diagnostic path, run existing focused tests and the full build, request review,
then run a fresh Task-16-equivalent B=6 calculation.
