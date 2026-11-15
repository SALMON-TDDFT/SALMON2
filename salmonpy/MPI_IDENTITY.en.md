# User IDs in MPI response sessions

This opt-in layer connects the serial point-run record to the existing two-level
MPI response service. It does not change integer transport IDs, rank allocation,
A/c-to-charge-current updates, reduction, or the field solver. Native backend
capabilities and restrictions still apply. There is no new MPI initialization on import.

All context ranks pass identical `point_record` and `point_map` to
`ResponseSession(..., dt=..., point_record=record, point_map=mapping)`, or to
`CoupledStepper.from_multiscale(..., point_record=record, point_map=mapping)`.
The record uses `salmonpy.point-run/1` (see [POINT_RUNS.en.md](POINT_RUNS.en.md));
it can be built as plain data or obtained from `PointRun.record()` without creating
serial states. Choose run_id on the host and share it before session construction;
do not independently generate UUIDs on every rank.

```python
record = {'schema': 'salmonpy.point-run/1', 'run_id': 'experiment-1', 'points': [
    {'point_id': 'right', **right_material.record()},
    {'point_id': 'left', **left_material.record()},
]}
# Existing plan.point_ids == (10, 90); the list and map explicitly agree.
mapping = {10: 'right', 90: 'left'}
session = ResponseSession(context, bindings, dt=0.02,
                          point_record=record, point_map=mapping)
```

Every integer point is mapped exactly once; ordered record entries must match
`plan.point_ids` through this map. Opaque user IDs are never sorted, normalized or
used as rank/file names. Material definition and optional initial-state hashes
are checked for internal consistency, and descriptors/records must agree across
all ranks before prepare/create callbacks. These are declared identities, not
independent proof that a provider actually loaded the named data: adapters still
must validate their source bytes and model compatibility. The cache key in
MaterialBinding is separate from material_id. A shared binding/cache key cannot
refer to different material or initial identities.

Providers and backends receive immutable `SolverContext.point_records`, keyed by
integer transport point for points in their own group. Each value includes run_id,
point_id, material_id/definition_digest and optional initial_state_id/digest.
Their existing point_id callback argument remains an integer. Without metadata,
point_records and `session.identity_record()` are None, preserving legacy usage.
With metadata, identity_record returns a detached mapping containing point_run and
mpi_point_map (transport point, user point, group name and world ranks). No new
state wrapper is inserted around collective native states.

## Collective output

Every context rank calls `session.save_arrays`. Only the field leader provides
arguments; others call it without arguments. Preflight and writer errors are
agreed across ranks. The field leader alone receives the OutputArchive.

```python
if context.group.role == 'field' and context.is_leader:
    archive = session.save_arrays(new_directory,
        {'current': current_history, 'time': time}, point_axes={'current': 1},
        metadata={'units': 'atomic', 'current': 'charge j'})
else:
    session.save_arrays()
```

Arrays must follow plan.point_ids order. Preserve timestamps explicitly; the writer
does not infer them. Reserved metadata keys: point_run, point_axes, mpi_point_map.
Read with `open_point_arrays`; it checks array hashes, point axes, user IDs and
transport/group ownership consistency. It does not restore MPI state or restart a
solver. All ranks must call operations in the same order. Process death and native
MPI abort remain unrecoverable. Use a bounded launcher/job timeout. A normal output
failure does not advance/invalidate the numerical session; retry with a fresh path.
The output directory is never overwritten. Metadata authenticity still requires
an externally trusted manifest hash.

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_point_identity /new/output
```

The sample also accepts 4 ranks (two ranks in one material group), compares against
a direct formula and the legacy session, and reads saved IDs back. It is synthetic
Python response, not TDDFT physics or a new FDTD solver. Use one thread per rank;
MPI/mpi4py must already be installed. No downloads or compilation are implicit.
