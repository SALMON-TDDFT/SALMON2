# R3: field initialization admission

The first recoverable stage is native E/H admission, before `read_input` and
numerical initialization. All field ranks participate even if one rank has a
zero point count, an empty/overlong path, or a missing input file. Rejection
must not consume `used`, create child communicators, or allocate field payload.
A corrected call on the same native attachment must initialize normally.

This does not make malformed namelist parsing, allocation failure, or failure
inside `eh_init` recoverable. Once admission succeeds the lifetime is consumed;
those later stages remain R3 work. The Python constructor's attachment/lifetime
policy is unchanged. The regression invokes rejected calls within that one
attachment before the normal constructor call, then checks normal evolution.

Validation: `tools/test_eh_deep_recovery_mpi.py --admission --inspect-mpi`,
1 and 2 field ranks, with rank-local invalid arguments. Compare all normal
A/E/H samples and output data bytes against the previous validated build.
No GS data is generated or copied. No GitHub Actions are scheduled.

Observed validation (2026-11-20): all three rank-local rejection cases followed
by valid initialization passed on 1 and 2 ranks; all 32 A/E/H samples and `.data`
outputs matched the R2 build exactly. MPI resource counts were zero after each
preflight rejection. Existing source tests: 638 run, 637 passed, 1 skipped.
Source/object/link prohibition audits passed. Record:
`calc2026/salmonpy-r3-field-admission-record-20261120` (canonical OneDrive results
root; cloud synchronization unconfirmed). R3 parser/deep failure work remains.

## Parsed admission (next stage)

Use an invocation-owned `s_input_parser`, so successful parsing closes its scratch
units before admission returns. A parsed request with non-Maxwell theory or
native restart enabled is rejected collectively before consuming the numerical
lifetime. The small context allocation and MPI-owner admission also precede
`used`; allocation failure is agreed collectively and allocated peers release
that empty context. `used` is set immediately before `eh_init`.

The regression repeats parsed restart rejection twice on the same attachment,
then initializes the ordinary input. This is native retry support; it does not
change Python constructor attachment ownership. Malformed parser input and
failures inside/after `eh_init` remain outside the recoverable scope.

Parsed-stage validation (2026-11-20): the previous build reproduced the consumed
lifetime failure. The updated build passed repeated parsed rejection followed by
normal initialization on 1/2 ranks. All 32 A/E/H samples and `.data` file bytes
matched R2 exactly, with zero child communicators after rejected calls. Source
tests: 638 run, 637 passed, 1 skipped. Source/object/link audits passed.
Canonical record: `calc2026/salmonpy-r3-field-parsed-record-20261120`.
No GS generated/copied; no Actions runs; cloud sync unconfirmed.

Distribution check: fresh sdist -> wheel -> clean environment succeeded;
installed tests also ran 638 (637 passed, 1 skipped). Installed solver source,
Japanese/English coupled guides and AI JSON exactly match the working source.

## Namelist read errors

`read_input`/`dump_input_common` accept an optional status used by native FDTD.
Detected positive namelist I/O statuses are broadcast and returned rather than
calling `end_parallel`/`stop`. The invocation-owned parser closes scratch files
on this rejection path. Callers omitting status retain standalone behavior.
Earlier unit selection, coordinate parsing, semantic validation and `eh_init`
failures remain separately scoped; this is not arbitrary parser exception safety.
The regression repeats an invalid integer in `&maxwell` then runs valid input.

Observed namelist-stage validation (2026-11-20): old build reproduced MPI
termination on the invalid integer. Updated build passed all seven admission
attempts (including two malformed integer attempts), then normal 32-step
propagation on 1/2 ranks. A/E/H samples and `.data` bytes matched R2 exactly.
Rejected calls retained zero MPI child communicators. Source tests: 638 run,
637 passed, 1 skipped. Source/object/link prohibition audits passed. Record:
`calc2026/salmonpy-r3-field-namelist-record-20261120`; no GS generated/copied,
no Actions runs, cloud synchronization unconfirmed.

Fresh sdist -> wheel -> clean-environment verification also passed: 638 tests,
637 passed, 1 skipped. Installed parser/field sources and Japanese/English guides
plus AI JSON match the working tree exactly. Changes remain uncommitted.

## Batch-job priority (user clarification)

Interactive input repair is not expected in a scheduled job. Prioritize checking
staged inputs on an allowed management/login node and recording rank-local
startup failures before solver creation. `preparation.preflight_inputs` reuses
existing integrity checks and an optional explicit trusted material registry.
CLI reports are exclusive-create files outside the immutable bundle; invalid
inputs return exit code2. `examples.mpi_input_preflight` checks each rank's local
view, gathers diagnostics, writes one report, broadcasts outcome and exits.
It does not edit inputs, submit jobs, execute solvers or retry automatically.
Runtime/native semantic admission remains separate, not certified by this check.

Verified valid management check, 1/2-rank checks, rank1-only missing input and
corrupted catalog rejection without solver startup. Record:
`calc2026/salmonpy-job-input-preflight-20261120`. No GS data collected and no
Actions runs; cloud sync unconfirmed.

Batch-input distribution verification: source and fresh wheel environments each
ran 640 tests (639 passed, 1 skipped). The installed preflight API, MPI sample,
Japanese/English preparation/API guides and AI JSON match the working tree.
Installed management check and MPI1/2-rank checks passed; rank1-only missing input
and corrupt catalog exited2 before solver startup. Record:
`calc2026/salmonpy-job-input-preflight-installed-20261120`.

## Ordinary mixed sample connection

`mpi_mixed_native --input-check STAGED_BUNDLE` reuses the MPI input-check helper
before native state creation and writes `OUTPUT/input-preflight.json`.
Optional `--input-manifest-sha256` anchors the manifest. This integrity-only gate
is opt-in and does not reinterpret the sample's legacy GS BUNDLE argument or
bind independently specified command arguments to staged material recipes.
Output-directory creation errors now broadcast before a barrier, preserving
existing outputs and ending all ranks with code2 instead of stranding peers.
Subsequent runtime/numerical failure policy is unchanged.

Source trial verified protected restore/continuation parity with/without input
checking (2 field ranks + 2 SBE ranks, LD local to field), and failed checksum/
existing-output collective exits before any solver creation. The initial test
harness expected `response.npz` from the stop-at-checkpoint driver; that driver
actually emits `protected-continuation.npz`. Existing outputs were retained and
compared using the correct artifact; no calculation code change was needed.

Fresh wheel normal sample: 4 ranks, 64 steps, A/E/j and point IDs exactly equal
with/without preflight. No GS collected/copied; no Actions run. Canonical records:
`calc2026/salmonpy-job-connection-20261120` and
`calc2026/salmonpy-job-connection-installed-20261120` (cloud sync unconfirmed).

User clarification: management/login-node preflight is an optional action the
user runs independently under facility policy. SALMONpy does not automatically
connect/run there or require a preflight report for ordinary job execution.
The mixed sample gate remains opt-in; no external scheduler integration is implied.

## Ordinary startup reads without any optional preflight

Moved explicit layers/model-group/TDDFT catalog reads before solver creation and
wrapped them in rank-local diagnostics, gather, report write and broadcast.
`startup-inputs.json` identifies rank, phase, path, exception type and message;
invalid declarations terminate all ranks with exit2. No management-node action,
report or --input-check is required. Detailed assembly/native semantic checks
remain later, with their existing policy.

Source and installed wheel verified missing layers, malformed JSON, rank1-only
read failure and parent MPI survival after caught rejection. Valid LD/SBE/native
FDTD protected continuation exactly matched the previous baseline. Both source
and installed tests: 640 run, 639 passed, 1 skipped. Installed sample and JP/EN/AI
guides exactly match source. No GS data collected; no Actions run; uncommitted.
Records: `calc2026/salmonpy-startup-reads-20261120` and
`calc2026/salmonpy-startup-reads-installed-20261120` (cloud sync unconfirmed).
