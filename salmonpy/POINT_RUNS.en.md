# Serial point identity and result arrays

`PointRun` binds an ordered mapping of point IDs to prepared models. It owns one
state lifetime per point, assigns an explicit run ID or generates a UUID, and
retains material and optional initial-state IDs/hashes. IDs are opaque,
case-sensitive Unicode strings, never used as filenames or inferred from ranks.

```python
from salmonpy.point_runs import PointRun, open_point_arrays

# model is a ConfiguredMaterial, including one returned by tddft_from_catalog.
with PointRun({'right': model, 'left': model}, run_id='experiment-1') as run:
    right = run.new_state('right', dt=0.02)
    left = run.new_state('left', dt=0.02)
    right.step([0., 0., 1.e-5])
    # Only right has advanced. LD takes E; SBE/TDDFT take A/c.
    identity = right.record()
    # Supply arrays ordered as run.point_ids == ('right', 'left').
    run.save_arrays(new_output, {'current': [[right.current, left.current]]},
                    point_axes={'current': 1},
                    metadata={'units': 'atomic units', 'note': 'points at different times'})
archive = open_point_arrays(new_output)
point_ids = [p['point_id'] for p in archive.metadata['point_run']['points']]
```

This wrapper does not schedule steps, align times, convert drives, infer units or
provide whole-run rollback. Record time arrays explicitly for physical analysis.
Each point can be created once per run, even after its state is closed. A failed
new_state call does not consume the point. Use a new run for branches or repeat
calculations. Run close attempts all state closes, re-raising the first failure;
failed closes are not retried automatically. Prepared models remain caller-owned
and must stay alive until all their states have closed. The wrapper expects
new_state to return a fully initialized state; it does not invoke its context-entry hook.
Resource-owning custom states must implement close(); context-only resource
managers are not supported. States without resources may omit close().
A state context closes just that point. Run close is idempotent. No concurrent
calls, MPI dispatch, checkpoints or restart mapping are provided in this layer.
Delegated snapshot/restore retains the backend's existing restrictions, without
adding identity to its native checkpoint format.

## User models

No inheritance is required. Supply `new_state(**kwargs)` and `record()` returning
exactly `material_id`, `definition_digest` (lowercase SHA256), plus optionally both
`initial_state_id` and `initial_state_digest`. Use a digest of the validated,
resolved definition, including relevant external-data hashes, not the object
address. Definition/initial identity is checked again before creating each state.
Keep other provenance in output metadata. Shared material IDs must have identical
definition digests; shared initial IDs must name the same material and initial
content. The same raw returned state object is rejected for multiple points, but
hidden shared arrays or native globals cannot be detected: independence remains
the custom backend's responsibility. No execution/physical capability is certified
by an identity record. Existing configuration registry validation still applies
when resolving JSON inputs; this Python binding is not an arbitrary JSON importer.

## Saving and loading

`run.save_arrays` adds an ordered `salmonpy.point-run/1` record under metadata
`point_run` inside the existing `salmonpy.outputs/1` archive. `point_axes` explicitly
maps array names to nonnegative point-axis indices. At least one array must have
a point axis of length `len(run.point_ids)`. Other arrays are left uninterpreted.
The caller must supply arrays in that exact order; matching lengths cannot detect
a permutation. The two metadata keys are reserved. Read back with
`open_point_arrays`: it verifies array hashes, shape/point mapping, IDs and identity
consistency. It never reorders data or restarts a solver. Manifest semantics are
validated, but metadata authenticity is not guaranteed; externally pin the manifest
hash if it will cross a trust boundary. Standard OutputArchive readers can still
read these files but do not validate the extra point schema.

```sh
python -m salmonpy.examples.point_ld /new/output/directory
```

This needs the installed LD wheel. It checks two independent points against the
direct LD path, saves currents/E/time and reads the point mapping back. It uses a
synthetic material, not experimentally validated parameters. No MPI or internet.

For explicit identity propagation/output in MPI response sessions see
[MPI_IDENTITY.en.md](MPI_IDENTITY.en.md). PointRun itself remains serial; restart is not provided.

Point-bound managed TDDFT disk restore is now available: [POINT_CHECKPOINTS.en.md](POINT_CHECKPOINTS.en.md).
The earlier restart limitations still apply to MPI whole sessions and general backends.

## In-memory snapshots with point identity

`PointState.snapshot()` wraps the backend-owned snapshot with immutable `salmonpy.point-snapshot/1` identity. `state.restore(saved)` checks run, point, material-definition and initial-state identity plus package version before calling the backend. Cross-point/run and raw backend tokens are rejected. Even with explicit matching run/point IDs, the backend must still validate its model, time step and other settings.

This is process-local rewind, not branching, cloning, disk persistence or an MPI-wide transaction. Unsupported custom backends raise `SnapshotError`. Backend restore-failure guarantees remain the backend's responsibility; identity rejection invokes no backend restore callback. Raw-state snapshot APIs are unchanged.

`python -m salmonpy.examples.point_ld /new/output/directory` also verifies cross-point rejection and exact current continuation after same-point restore.
