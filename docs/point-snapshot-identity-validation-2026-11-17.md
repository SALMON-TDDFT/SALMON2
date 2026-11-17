# Point snapshot identity — 2026-11-17

## Scope and regression

PointState previously forwarded snapshot/restore directly to its backend. A custom backend accepting an opaque tuple reproduced cross-point restoration without any point identity check. The failing regression was observed before implementation.

The wrapper now returns an immutable process-local point envelope containing exact run/point/material-definition/initial-state identity, schema and package version. Restore checks these before invoking the backend. Opaque payloads pass unchanged to their originating backend protocol; backend model/settings checks and restore-failure guarantees remain its responsibility. This is not a security boundary, disk format, cloning/branch API or MPI transaction.

## Verification

- Source suite, with explicitly selected built LD/SBE library: **517 tests passed, one skipped**.
- Built sdist then wheel, installed in the dedicated Python3.14 Apple Silicon environment: **517 tests passed, one skipped**.
- Point tests cover a custom opaque tuple payload, callback-count rejection before mutation, every identity field, schema/version, immutability/detached records, closed/unsupported states and raw-token rejection.
- Actual LD/SBE states: three future steps exactly match after same-point rewind; cross-point rejection leaves the other point unchanged. An additional targeted run of all10 point tests confirms matching identities still reject a different dt through backend validation, preserving its state.
- Installed `python -I -m salmonpy.examples.point_ld NEW_DIRECTORY` passes direct reference, independence, cross-point snapshot rejection, exact restored continuation and identity output roundtrip. Original20-step output endpoint is retained.
- Read-only code review found no blocker. JP/EN guides and AI JSON copies match; whitespace check passes.

No numerical formulas, solver binaries, native payloads or MPI orchestration were changed. Only local tests/builds were used; no GitHub Actions consumption. Dedicated evidence directory: `/Users/otobetoshihito/.codex/outputs/salmonpy-point-snapshot-identity-20261117/`. Local paths record provenance and are not user requirements.

## Remaining

Other backend disk formats, explicit branch identity, changed-MPI-topology redistribution, unrestricted native context ownership and recoverable native errors remain open. Memory estimates remain warnings only. Molecular work remains paused.
