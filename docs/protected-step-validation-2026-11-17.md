# Common protected step — 2026-11-17

## Scope

Added public preview StepTarget, guarded_step and StepRecoveryError. Each target
binds an owner to explicit capture/rollback/release callbacks. Backups stay opaque
and local; no copying, serialization, MPI token transfer, native capability
inference or propagation replay. Solver-group key order and optional release
availability agree before capture. Global synchronized slots support uneven
numbers of targets. Partial capture cleanup calls release(None) on the failed
rank, preserving the group's collective schedule. Rollback and cleanup continue
after an individual callback failure; preservation is never claimed if any
rollback fails. Cleanup failure after successful operation may leave a step
accepted and does not claim original preservation.

CoupledStepper.step_protected requires exactly one explicit target per local
field/material object and protects coupling clocks/current/A/adapter/measurement
metadata separately. Incomplete rollback/cleanup poisons the response and local
field. Ordinary step remains unchanged. Contracts must include all backend
physical state and external side effects; declared methods are not certification.

This first implementation covers explicitly registered light/custom components.
It does NOT supply aggregate native E/H/TDDFT/fast-PR HSE numerical/output
transactions. Their individual recovery protocols are insufficient for this
purpose. Native outputs and handle/lifetime side effects need a separate adapter;
no native targets or unsafe fallback are generated automatically. Process/MPI
termination, native STOP, deadlock and allocation exhaustion remain excluded.

## Verified locally

Apple Silicon, user host64GB/18 cores, six ranks with OMP/BLAS1. No remote Actions,
solver source/binary changes, automatic runtime build or network dependency fetch.

- Six-rank example joins Multiscale1D1 + real LD1 + real SBE1 + a three-rank custom
  group with two independent macro points. Target counts differ across groups.
- Late one-rank custom mutations/failures at accepted attempts1/25/80, LD mutation
  at40 and failure after field mutation at50: all numerical snapshot bytes,
  wrapper A/current/clocks and coupling histories/flags are unchanged afterwards.
  Original coupling A/J/adapter arrays and example-wrapper A/current references
  are preserved too, not replaced by equivalent backup arrays.
  Valid explicit retry produces80 accepted frames exactly equal to ordinary
  execution; zero elapsed propagation replay. The custom fault trigger is test
  control, deliberately outside the physical-state token.
- Missing owner rejects before capture/advance. Mismatched optional collective
  release rejects before ANY registered capture (audit list remains empty).
  One custom rank's partial capture failure still calls collective release(None),
  preserves originals and permits retry. Rollback failure and success-followed-by-
  release failure invalidate every group's response and report unsafe status.
- The installed-wheel six-rank module example passes the same ten cases without
  source path injection. Source and installed suites:534 tests passed, one skip
  each. Six dedicated serial tests include opaque/None token ownership, partial
  capture cleanup, empty-target admission and unsafe rollback/release reporting.
- Source distribution -> wheel, native hashes, packaging layout, installed manual
  search and Markdown/AI JSON synchronization pass. Existing dependencies were
  reused offline from the local test environment.

## Review and negative evidence

Initial unit test failed before stepping.py existed. Review found optional
release availability was not initially part of subgroup admission; it is now
agreed with key order and covered by a zero-capture real MPI negative. Final
read-only review found no remaining blocker within the explicit contract. The
alias-strengthened regression also exposed metadata copies replacing accepted
array references; originals now restore in place and the strict identity checks
pass. The first full suite caught missing environment/source fields in the new AI catalog
task; those metadata were added and full source/installed suites now pass. These
fail-first records are retained, not counted as successful verification.

## Reproduction

`mpiexec -n 6 python -m salmonpy.examples.mpi_protected_step NEW_OUTPUT`, with
OMP_NUM_THREADS=1, BLAS1 and an external timeout. This is a synthetic implementation
regression, not physical material convergence. All outputs must use a new directory.

Persistent local artifacts:
`/Users/otobetoshihito/.codex/outputs/salmonpy-common-protected-step-20261117`.
Includes wheel/sdist, source hashes, logs, installed results and accepted arrays.
