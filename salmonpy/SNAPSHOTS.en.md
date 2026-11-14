# Saving and restoring lightweight states

[日本語](SNAPSHOTS.md) · Topic: `snapshots` · Paired revision: `2026-11-14.1`

`state.snapshot()` returns an immutable `StateSnapshot` token.
`state.restore(token)` replaces the accepted state and returns `None`.
These operations are process-local, fixed-dt, in-memory operations. They do not
execute a numerical step, read/write a file, invoke MPI or compile/load a solver.

## Supported history

| State | Saved evolution history |
|---|---|
| `LDState` | Per-pole P/J, previous E and step counter |
| `SBEState` | Density matrix, previous A/c, cached native matter current and counter |
| `Transverse` state | A and step counter |
| `Longitudinal` state | A, two induced-A levels and counter |
| `ThinFilm` state | A, reflected A and counter |
| `Multiscale1D` state | Macro A, full old/current/new grid including ghost cells, grid current and counter |

The token also records schema, exact package version, time, compatibility and
payload integrity. It retains the originating model object; discard unneeded
tokens to release memory. Snapshots copy evolving arrays, so memory use grows
with saved state size and snapshot count. There is no compression or memory cap.

## Compatibility and isolation

Restore accepts the same state or another state created from the **same model
object**, with identical dt and time origin. LD also requires identical initial
field configuration. Fields require identical capacity, geometry/boundaries and
sampled external/incident drive. An equivalent but separately constructed model
is deliberately rejected in this first version. For a callable drive, each new
field state evaluates the drive; differing samples are rejected even when the
callable/model identity matches.

Validation and array allocation finish before mutation. Invalid tokens raise
`SnapshotError` and preserve all accepted history. Tokens contain immutable byte
copies; restoring creates new arrays. Advancing/restoring one branch does not
change another branch or an existing token. Built-in state objects themselves
are not safe for simultaneous calls from multiple threads.

```python
import numpy as np
from salmonpy.ld import LorentzDrude

model = LorentzDrude(plasma_frequency=1., strength=[1.],
                     resonance=[0.], damping=[0.1])
a = model.new_state(dt=0.02)
a.step([0.1, 0., 0.])
saved = a.snapshot()
expected_p, expected_j = a.step([0.2, 0., 0.])
a.restore(saved)
p, j = a.step([0.2, 0., 0.])
np.testing.assert_array_equal(j, expected_j)

b = model.new_state(dt=0.02)
b.restore(saved)             # independent branch at saved.time
b.step([-0.1, 0., 0.])      # does not change a or saved
```

Run `python -m salmonpy.examples.snapshot_restore` for LD, SBE and all four
fields. It needs installed LD/SBE kernels but no full SALMON runtime. The sample
checks exact replay locally; it does not claim cross-platform bitwise identity
or physical convergence. Ordinary scientific regression tolerances remain valid.

## Explicit exclusions

This is not the handoff's Zarr GS archive, an output NPZ, a persistent restart,
a cross-version/cross-process format, a C ABI or a variable-dt method. Do not
serialize/edit tokens or pass output archives to `restore`. Snapshot public
metadata (`kind`, `schema`, `version`, `step_index`, `time`) is diagnostic only;
private fields are not an editable checkpoint schema.

TDDFT workers, embedded TDDFT, native E/H FDTD and MPI sessions do not gain this
API. A host manually coupling supported material/field states must save and
restore **both at the same accepted coupling boundary**. There is no automatic
transaction across components; a single state's atomic restore does not provide
MPI/global rollback. `couple()` still creates and closes its own states.
Changing dt, adding adaptive stepping/subcycles or making native TDDFT reentrant
requires separate work. Transverse prescribed-A use remains the basic workflow.
