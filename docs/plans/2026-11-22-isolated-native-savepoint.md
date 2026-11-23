# Isolated worker native savepoint implementation plan

Goal: implement the first layer of no-replay molecular rollback in the existing
PZ E worker, preserving the public E3 numerical reply and normal calculations.
This phase is a live-worker savepoint, not a durable disk checkpoint or Python
public snapshot replacement. The existing replay API remains explicitly replay.

Architecture: caller-owned savepoint local to main_tddft (one pending backup),
containing both orbital buffers, temporary orbital, densities/potentials,
Poisson state, RT histories, energy, field vectors and XC payload. Explicit
CAPTURE_E/RESTORE_E/RESUME_E/RELEASE_E commands; output streams close before Python
journal truncation and reopen after. Status replies include accepted step and a
monotonic count of actual time_evolution_step calls, never restored. No global
selector/cache, MPI initialization or extra physical reduction. Unsupported
output modes reject only savepoint creation, not the underlying calculation.

1. Extend the existing raw isolated-worker validation with capture/continue/
   restore/resume, counter and process identity assertions, repeated restore,
   invalid-state rejection, independent worker and unchanged output checks.
   First verify failure against the unchanged E3 worker.
2. Implement worker command dispatch and a main-local numerical savepoint. Use
   an explicit loop clock so restoration can select the saved parity/index.
   Ordinary batch/periodic stepping order and numerical calls remain unchanged.
3. Update byte-exact inventories; prepare/build a tracked-only source export,
   retaining source/object eigensolver checks. No new GS; reuse the authorized
   C2H2 campaign. Test initial, odd/even and capacity endpoints, all observations,
   sparse energy freshness, stream rollback and no extra evolution calls.
4. Run NLCC local MPI1/2 tests before advancing the source capability identity;
   run source/wheel suites and existing public coupling regression.
5. Document evidence and remaining public Python/disk integration separately.
   No claim of arbitrary allocator/process failures, other XC, MPI-internal
   isolated states or portable persistence in this first phase.
