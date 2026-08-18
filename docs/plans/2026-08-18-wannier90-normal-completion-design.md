# Wannier90 Normal Completion Design

## Context

Wannier90 can reach `num_iter`, write its final state and checkpoint, print its
normal completion banner, and return successfully without printing the optional
early-convergence message.  SALMON currently rejects that run before inspecting
the returned transformation because its convergence-log parser requires the
early-convergence marker.

## Decision

Treat a readable `.wout` containing both a final state and Wannier90's normal
completion banner as a completed Wannier90 run.  Preserve an early convergence
iteration as a diagnostic when present; otherwise report `num_iter`.  Continue
to reject missing, truncated, or abnormal logs.

After this log-level completion check, retain all existing SALMON result gates:
finite transformation/centres/spreads, unitarity, spread checks, and downstream
translation and point-cogroup symmetry validation.  Thus iteration exhaustion
alone is no longer fatal, while malformed numerical output remains fatal.

## Testing

The MPI fixture will prove that a normal final state at the configured iteration
limit is accepted and reports that limit.  It will also retain the early-
convergence case and add an abnormal/truncated log rejection so the change does
not weaken failure detection indiscriminately.
