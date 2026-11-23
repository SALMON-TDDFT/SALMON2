# Atomic restore through isolated candidates

User scope: preserve original state on mid-restore failure, including MPI/custom backends; update paired/AI manuals.

- Add CheckpointState, one process-local holder for field and all owned point states. A single internal bundle pointer publishes all local candidates together; old references stay alive through commit. Read field/states through the holder after a successful restore.
- Extend CheckpointTarget with stage_restore() -> independent fully loaded state and optional discard_restore(candidate_or_None). Staging validates/loads privately and never mutates live states. Partial allocations on a throwing stage callback remain that callback's responsibility.
- Add restore_checkpoint as the protected entry, and atomic=True/live= on execute_checkpoint. Existing execute_checkpoint default remains non-atomic for compatibility. Atomic mode rejects unsupported direct-only targets before any callback; no silent fallback. Built-in in-process native contexts are not cloneable through the current ABI and remain explicitly outside protected mode.
- All ranks admit, reserve holders, stage slots, agree failures and vote before one non-callback bundle assignment. Reject reused live/candidate identities. Abort discards attempted candidates in reverse synchronized slots, including None for a stage which threw. Cleanup failures are reported separately; live holders remain unchanged.
- AtomicRestoreError(RuntimeFailure) reports original_preserved and cleanup_errors. Normal exceptions before publication preserve the live bundle; fatal process/MPI failures and backend violation of candidate isolation cannot be recovered generically.
- Verify multi-component later-stage failure, commit success, alias/unsupported admission, cleanup failure, continued original evolution, real MPI unequal groups/one-rank failure and custom JSON plus LD/SBE candidate restore. Recheck existing native non-atomic paths unchanged and packaging/guides.
