# Point-bound TDDFT checkpoint validation (2026-11-15)

M6c: PointState disk checkpoint methods wrap the existing managed TDDFT replay
archive in point-checkpoint/1. The outer envelope contains exact run/point/material/
initial IDs and hashes, version, and SHA256 of the inner manifest. Identity mismatch
is rejected before native replay. No model/runtime is loaded from archive paths.
No changes to TDDFT numerical propagation, internal archive layout or replay logic.

Eight focused tests cover same-identity roundtrip, different point/run rejection
before restore, manifest corruption, extra entries, duplicate/oversized JSON,
failed-save cleanup/no overwrite, unsupported backend, invalid byte budget,
symlinks and closed states. Review found no important static-input defect and
prompted the extra bounds/lifetime/path tests and explicit concurrent-modification
restriction. Native input and outer manifest are separate reads, not an authenticated
race-resistant transaction.

Source suite: 389 tests, one filesystem-condition skip.
Fresh installed wheel suite: 389 tests, two skips (filesystem and optional plot).
sdist -> wheel packaging passed; source-independent imports used with python -I.

## Real installed-wheel Si replay

A fresh packaged PZ Si GS was generated and transferred to a pinned GS/1 catalog.
Sample examples/point_tddft_restart.py created point p, advanced 10 steps, saved
its checkpoint, and continued to step 20. Restore into point other was rejected
without changing its t=0 current/time. The first run and prepared model were closed.
A newly loaded model with the same declared initial data and a new PointRun with
the same explicit run_id/point_id restored the checkpoint and continued 10 steps.
All resumed currents matched the uninterrupted history exactly (maximum absolute
difference 0). IDs remained Si-starter / Si-ground / p / resume-example.

Output/report: `/tmp/salmonpy-point-checkpoint-native-20261115`.
Input preparation/report: `/tmp/salmonpy-point-checkpoint-inputs-20261115`.
Logs: `/tmp/salmonpy-point-checkpoint-*.log`.
Generated GS restart payloads were removed after verification to avoid accumulated
copies; report, checkpoint history, manifests and current arrays remain. Regenerate
GS inputs to rerun the example. Checkpoint does not itself contain GS wavefunctions.

macOS Apple Silicon, Python 3.14; existing managed runtime executable SHA256
29444c0e75447a63611d90b2811207a2698c76a4a18005169543cc6c2f039bed.
One thread/worker, maximum two simultaneous workers. Peak memory not measured or
hard-enforced. No MPI run, native build, GitHub Actions, commit or push.

## Limits

This is short PZ Si connection/replay evidence, not optical convergence or other-XC/
cross-platform evidence. Exact IDs are required for resume; cross-run branching
is deferred. LD/SBE/custom disk formats, MPI field/material/global restart and
in-memory point identity wrapping remain separate work. Direct native checkpoint
APIs and in-memory snapshots keep their original contracts.

Envelope JSON is bounded at 64 KiB, inner manifest prehash at 16 MiB (and any smaller
positive max_bytes). Native loader retains its format-derived limits and atomic
replay. max_bytes is an inner archive/payload budget, not RSS. Keep archives fixed
while reading and authenticate outer metadata separately when crossing trust boundaries.
