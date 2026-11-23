# Serial point/run identity validation (2026-11-15)

Implemented the serial subset of approved M5/M6 in salmonpy.point_runs:
PointRun, PointState and open_point_arrays. Prepared models remain borrowed.
PointRun owns at most one state lifetime per point, records exact opaque IDs,
and saves an ordered point mapping plus explicit point axes in outputs/1 metadata.
No numerical kernels, field scheduling, current signs or native input are changed.

Seven focused tests cover two-state isolation with a custom Python backend,
identity retention and detached records, invalid IDs and conflicting definitions,
failed-creation retry, exact shared-object rejection, cleanup of all states,
borrowed models, changed identity, reserved metadata, point axis validation and
archive roundtrip/malformed point mapping rejection. Review clarified that custom
new_state must return an initialized state and resource owners must implement close;
context-only resource management is unsupported. Hidden array/native-global sharing
cannot be certified by this wrapper.

Local source suite: 367 tests, one filesystem-condition skip.
Isolated installed wheel suite: 367 tests, two skips (filesystem condition and
optional plotting dependency). Python 3.14 / Apple Silicon / macOS 26 only.
Build uses sdist -> wheel, at most two jobs, and explicit local dependencies.
No GitHub Actions or full solver runtime rebuild was performed.

The native example `python -I -m salmonpy.examples.point_ld NEW_DIRECTORY` ran from
an installed wheel outside the source checkout. Two points of the same synthetic
LD material receive opposite E drives for 20 steps. Updating right leaves left
unchanged; right agrees exactly with an independent direct LD reference at every
step. Saved current arrays and material/point/run identity roundtrip exactly.
The synthetic model is not a validated physical material. The test does not prove
MPI independence, custom-backend hidden-storage isolation or physical convergence.

Reports/arrays: `/tmp/salmonpy-point-runs-installed-example-20261115`.
Build/test logs: `/tmp/salmonpy-point-runs-*.log`.
No hard memory-budget enforcement or peak-memory claim. Builds limited to two jobs;
native example uses one process and one thread. Existing sources and runtime stay
unchanged. No commit/push requested in this continuation.

Remaining: MPI ID/placement integration, state-restoring checkpoints/restart ID
compatibility, explicit backend capabilities, full native state isolation. Output
mapping validates declared order/shape, not whether caller supplied the correct
physical permutation or trustworthy metadata; pin manifest separately as needed.
