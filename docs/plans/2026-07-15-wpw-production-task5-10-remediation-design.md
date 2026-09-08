# WPW Production Task 5-10 Remediation Design

## Scope

Complete Tasks 5-10 of the accepted Wannier+PW DG-DC and Exp plan.  The
checkpoint ends after a Si64 2x2x2-fragment DG-DC, field-off RT, and 20-step
laser smoke test.  Tasks 11-12 and long HHG production runs remain out of
scope until the checkpoint review is complete.

## Production data flow

The production basis is distributed `windowed_kg` with stable column id
`(K-1)*n_G+G_id`.  Rank-local volume and canonical-face quadrature produce
owned sparse WP/PP candidates.  A provider/halo layer supplies remote face
traces; the canonical scanner performs no communication and publishes a face
block only after a complete successful scan.  Sparse blocks feed matrix-free
H/S callbacks and never form global dense H, S, density, WP, or PP matrices.

DG-DC retains only the occupied subspace plus a bounded extra-state window.
It reconstructs the complete mixed density on uniquely owned core points,
uses the existing distributed Hartree and LDA paths, and converges density,
potential, energy, occupied projector, generalized residual, metric
orthonormality, charge, and gap gates.  The existing dense fixed-basis solver
remains a small mathematical oracle and is never called by `main_dft`.

## Checkpoint and RT handoff

One GS/RT-neutral checkpoint schema stores basis/operator fingerprints,
distributed layout, metric policy, retained eigensystem, occupations,
converged potential, fixed operator blocks, face convention, and checksums.
RT loads this state without projecting conventional orbitals or reselecting
the metric.  Before propagation it verifies rank-local H/S block identity and
occupied generalized residuals.

## Propagation and observable

The accepted production propagator is a namelist-controlled midpoint Exp in
the S metric.  Each midpoint correction restarts from saved `C_n`; failure to
converge or preserve the S norm is terminal.  Field-off stationarity is a hard
gate before field-on execution.  Length-gauge WPW blocks follow the accepted
operator contract, retain periodic-H1 PP with no PP face term, report the raw
metric-Hermiticity residual, and track a continuous `Delta_Pz` branch with
`Jz=d Delta_Pz/dt` as a secondary consistency signal.

## Si64 checkpoint

Task 10 creates a provenance manifest and matched Si64 2x2x2-fragment DG-DC
and RT inputs.  The execution order is DG-DC convergence, checkpoint identity,
field-off RT, then a 20-step laser smoke.  No long HHG spectrum is interpreted
at this checkpoint.  After Task 10, review Tasks 5-10 findings-first and add a
RED regression before every review-driven production fix.

## Safety and validation

Every production change follows RED-GREEN-REFACTOR and is committed as a
small task checkpoint unless the user explicitly requests no commit.  Failures
return no partial blocks or checkpoints.  No hidden MPI, global fragment scan,
O(N) owner metadata, environment-only scientific control, dense production
operator, push, or pull request is permitted.
