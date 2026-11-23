# Python live savepoints and molecular cancellation

Goal: expose the verified worker savepoint and use it for FDTD cancellation with
no electron replay. One pending live-worker savepoint, explicit lifetime; replay
archives/disk checkpoints remain explicitly separate. No compatibility fallback.

Add IsolatedSavepoint and state.savepoint/restore_savepoint/release_savepoint.
Bind tokens to the actual engine and worker PID, preserve immutable observation
metadata, truncate accepted Python histories after native restore and retain the
worker's monotonic update count. Use the existing output journal between native
RESTORE_E and RESUME_E. A journal failure suspends stepping and allows explicit
restore completion; no automatic process replacement. When a live savepoint is
active, step failure never starts a replay worker. Unknown/dead workers are
unusable and whole-coupling recovery must report original_preserved=False.

Use current validation receipt revision4 requiring native savepoint proof and
same-process/no-extra-evolution evidence. Reject revision3, with no version rescue.
E3 numerical reply stays nine values; add the state reply prefix to the existing
reader queue. Current replay archive format is unchanged, not inferred.

1. Extend the existing synthetic worker/receipt and add a failing lifecycle test:
   same PID, zero restore updates, independent state, stale/foreign token rejection,
   busy replay rejection, repeated direct restore, explicit release.
2. Implement in isolated_tddft.py; integrate release on accepted/cancelled molecular
   trials. Keep failed restoration tokens until restoration can be completed.
3. Validate real C2H2 public savepoints against old continuous trajectories, then
   MPI6 FDTD reject30/retry: same electron PID, one rejected evolution only,
   restore adds zero updates, all accepted arrays and six-face budgets exact.
4. Run source and installed-wheel tests/examples; update paired guides, AI JSON,
   checklist and reproducibility evidence. No new GS/optional physical CI.
5. Durable direct disk persistence and process-loss restoration remain separate.
