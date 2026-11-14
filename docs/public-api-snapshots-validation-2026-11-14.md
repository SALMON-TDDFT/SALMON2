# Public API and state snapshots validation

Date: 2026-11-14. Local macOS arm64, Python 3.14. Package 0.1.0a5.
This implements the user-approved lightweight, fixed-dt, in-memory checkpoint
scope after reading docs/plans/2026-11-14-handoff.md. Full solver source import,
Zarr GS archives, C ABI and native RT reentrancy remain separate work.

## Behavior verified

- Public version/inventory, quiet import with optional integrations unloaded,
  distribution-version agreement, and legacy built-in exception catches.
- LD replay with multiple poles, damping and nonzero initial E.
- SBE replay with all four correction orders and both nonlocal-correction modes;
  nonzero time origin/initial A and original native-reference tests retained.
- All four fields, plus ABC incidence, PEC and periodic Multiscale1D boundaries.
  Full grid A/E as well as macro state are compared after restore.
- Restore into a separate compatible state, no shared evolving arrays, immutable
  saved bytes, initial-state and exhausted-field rewind, altered-drive/model/dt/
  time-origin/capacity rejection, corrupt payload/metadata rejection and unchanged
  state after failure. A manual SBE/Longitudinal coupled loop also replays exactly.
- These are local exact-replay checks; no cross-platform bitwise or physical
  convergence guarantee is inferred. Original numerical update bodies are retained.
- Python source suite: 171 tests passed. Installed-wheel suite: 171 tests passed.
  The installed test runner rejects checkout imports; the import subprocess also
  verifies it loaded the same package. The installed six-component snapshot sample
  passed, and AI manual search works from the wheel.
- Source distribution -> wheel build passed; installed native hashes passed.
  New Python modules and paired API/snapshot guides are present in the wheel.
  `tools/sync_manual_guides.py --check` passes; every bundled guide is compared
  with its embedded JSON text by a test.
- Independent code review found no material correctness issues.

No changes to native/ or runtime/source.json/patches were needed. No full-native
MPI runtime regression or new platform CI run is claimed for this work. Test
runs are small single-process jobs; no new peak-memory or hard-limit measurement
was made. Snapshot storage copies state and retains the model; users must budget
memory for each retained token.

Validated wheel SHA256: `968877cc187b7b8cf86c341727d9499fca173cd15198f6c60b8bb46a4eda965c`.

## Distribution contract

The public API remains preview. Explicit package errors use SalmonpyError with
input/numerical/runtime families; foreign OS/dependency/callback/conversion
errors may keep their original types. RuntimeUnavailable/OutputSaveError retain
their names. Snapshot restoration requires the exact package version and same
model object, compatible fixed dt/time origin and field drive/capacity. There is
no durable serialization, adaptive dt, TDDFT/native-EH/MPI restore or atomic
multi-component transaction.
