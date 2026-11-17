# Protected checkpoint restore — 2026-11-17

`restore_checkpoint` stages independent backend-owned candidates and publishes a
single `CheckpointState` bundle only after every MPI participant reports success.
A handled staging failure discards attempted candidates in reverse slot order,
retains the original bundle, and raises `AtomicRestoreError`. Cleanup errors are
reported separately without replacing the original bundle. Successful calls
return previous objects for explicit retirement by the host.

## Evidence

- Source and installed-wheel suites: 495 tests pass, one skipped in each run.
  Actual LD and SBE candidate restores reproduce future currents exactly. Tests
  cover late failure, original-object rejection, direct-only rejection, cleanup
  failure, independent states and old-state continuation.
- Real four-rank MPI: field rank, one-point rank, and two-rank/two-point material
  group. A rank fails after loading and modifying its private candidate at the
  later point. All original references, arrays, clocks and counters remain
  unchanged; old states continue correctly. Retrying successfully loads step17
  and advances15 steps to the continuous step32 result exactly. Inconsistent
  atomic policy/discard availability and direct-only callbacks are rejected
  before loading. The host communicator remains usable.
- The same MPI test passes against the installed wheel. Existing four-rank
  checkpoint execution and three-rank restart-preflight tests also pass.
- Installed custom JSON example reports `failed_restore_original_preserved`,
  independent states, zero maximum continuation error and model-mismatch rejection
  before candidate loading.
- Source-distribution-to-wheel build and installation pass. Japanese, English
  and AI JSON guides describe the same ownership and failure contract.

Persistent logs and built artifacts:
`/Users/otobetoshihito/.codex/outputs/salmonpy-atomic-restore-20261117`.
Tests used at most four MPI ranks with OpenMP/BLAS threads set to one. No remote
Actions were run and no native kernel or numerical update formula was changed.

## Limits and remaining work

Factories must not mutate or share mutable storage with live objects. A failed
factory owns its partial allocations; `discard_restore(None)` must be safe and
complete any required subgroup communication. The API cannot rescue stranded
collectives, process crashes or hidden backend aliasing. Hosts pause evolution
during restoration and retrieve new objects from the holder after success.

The existing native in-process FDTD/TDDFT/HSE ABI cannot construct isolated
candidates. Protected mode rejects these direct-only targets before LOAD;
the supplied native coordinator explicitly retains non-atomic execution. Native
atomic candidate support remains open. This validation does not claim native
rollback, changed-topology redistribution or cross-platform numerical validation.
