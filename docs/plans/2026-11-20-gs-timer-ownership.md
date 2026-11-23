# Group5: GS-owned rank-local timer diagnostics

Split GS timing storage from process-wide timer arrays. Each retained GS owns a
parked timer context; create/continue borrow that context, restore the host bank
on all normal returned paths and exclude parked intervals from running timers.
Move allocations instead of copying numerical arrays. Do not add timer data to
SCF numerical payloads or change arithmetic, mixing or native stream1.

Expose optional native timer protocol1 and rank-local completed interval seconds
through NativeGSState.timing_seconds. Old runtimes retain normal GS operations;
only the new query requires the additive capability. Timers are lifetime-local
diagnostics, not elapsed-time replay or wall-time comparisons across MPI ranks.
Fresh-process restart starts fresh timers; numerical restore does not rewind
already spent work. Threaded concurrent native execution remains unsupported.

Verify deterministic parked/running/thread/subrange/host/sibling cases with
fixed-clock native fixture, including old-source compile rejection; verify live
ordinary/full/reduced HSE owner isolation, protected restore and previous-runtime
bit equality using the existing MPI test (numerical arrays in RAM, no GS export).
Run source/wheel tests and package the native fixture/runner, human JP/EN and AI
instructions. No routine CI dispatch. Keep scalar evidence in calc2026.

Later group5 work remains shared I/O/initialization/MPI ownership, wider stable
ABI coverage and recoverable allocation/process/publication failures. This does
not establish concurrent or arbitrary fatal-error recovery.

## Verification

Deterministic source and installed native fixtures pass. Ordinary MPI2 and
full/reduced-k HSE MPI4 with three independent groups: previous-runtime physical
arrays bit exact, sibling timer dictionaries unchanged, protected failure and
retry exact, restore SCF work0. Wheel ordinary/reduced-HSE controls pass.
Source and installed tests637:636pass/1skip. Numerical stream writer bytes are
unchanged; diagnostics are not archived numerical state. Source/object audits
passed. No GS datasets collected and no CI dispatched. Scalar evidence lives in
calc2026/salmonpy-gs-timers-record-20261120; cloud sync is unconfirmed.
