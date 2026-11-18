# Context-owned RT transactions and restore candidates

Next step toward stable native handles: remove persistent process-global protected
step buffers, pending restore candidates and their output stream metadata from
salmonpy_tddft_context. These must belong to the actual tddft_context. Concurrent
transactions on one context still conflict; independent contexts must not block,
replace or discard each other's data. Close deallocates only its own candidates.

First verify the existing internal non-symmetric fast-PR route can open two contexts
sharing immutable physical settings and the same MPI communicator. Compare isolated
and interleaved calculations with different dt/step capacity and waveforms. Then
hold simultaneous step transactions, restore one while the other advances, stage
both restore candidates, commit/cancel separately, and close one with the other's
candidate still alive. Repeat exact/learned Taylor/PT-CN on internal MPI1/2.

No numerical formulas or checkpoint encoding changes. Keep existing scoped output
close/reopen and host file rollback contract; snapshots do not intrinsically copy
owning FFT/C/MPI handles. Module-scope serializer allocation bound is synchronous
scratch, not a pending candidate; leave it explicit as a remaining extraction.

Public one-lifetime/backend and immutable common input/symmetry/MPI requirements
remain unchanged. This is preparation for a generation-checked C handle registry,
not a claim of arbitrary heterogeneous same-process execution. Native communicator,
I/O/timer/remaining common input ownership and handle ABI follow after verification.
Run source/installed suites, own restart/mixed regressions, JP/EN/AI synchronization
and sdist-to-wheel after native layout changes. No Actions.

Review also identified closed output-unit reuse: opening B during A's rollback
window can take A's former positive unit numbers. Reopen must allocate available
units and update A's ofile; closing an already rolled-back A must not close B.
Add late-open/reopen and late-open/close cases and a small direct output probe.
