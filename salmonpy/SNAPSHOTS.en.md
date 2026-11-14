# Saving and restoring states

[日本語](SNAPSHOTS.md) · Topic: `snapshots` · Paired revision: `2026-11-14.3`

`state.snapshot()` returns an immutable `StateSnapshot` token.
`state.restore(token)` replaces the accepted state and returns `None`.
These operations are process-local, fixed-dt, in-memory operations. They do not
advance the accepted time. Lightweight restore copies arrays; process-isolated
TDDFT restore starts a fresh audited solver and replays accepted steps using
private temporary files. Neither API builds a solver or invokes MPI.

## Supported history

| State | Saved evolution history |
|---|---|
| `TDDFTState` (process isolated) | Initial A/current, accepted A/current history and counter |
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

## Process-isolated TDDFT replay

Use `saved = state.snapshot()` and `state.restore(saved)` on a state returned by
`TDDFTInitialState.new_state`. The same initial model object, fixed dt, initial A,
capacity, executable checksum and thread count are required. Closed states reject
both operations. Existing states retain a lease on one shared baseline after
the model loader closes; each worker owns an independent writable copy. N live
states retain N+1 managed bundles, excluding caller input and saved outputs.
Restore temporarily adds one worker copy, sharing the same baseline. Closing
the loader refuses new states; releasing the final state then removes the
baseline. Every spawn verifies the copied complete manifest before process
launch. Worker file changes do not propagate to the baseline or other workers.

Restore starts a separate worker and replays every accepted A, checking initial
and endpoint matter currents with rtol=1e-12, atol=1e-14. Only after validation
succeeds is the old worker replaced. A replay failure raises `SnapshotError` and
preserves the original live branch. The actual replayed current is retained.
The accepted branch's output files are regenerated, so a later `finish()` exports
the restored branch. Discarded trial outputs are not merged into it.

Snapshot costs O(saved steps) memory; restore costs a fresh initialization plus
O(saved steps) evolution and temporarily holds two workers and their files.
The token does not store wavefunctions and cannot survive a Python restart.
Same-machine single-thread tests check bitwise continuation; cross-platform or
parallel bitwise equality is not promised. Coupled hosts must also rewind fields.

## TDDFT disk checkpoints

```python
# state is an existing TDDFTState. The destination must be new.
state.save_checkpoint('checkpoint-0005')

# In another Python process, explicitly supply the identical initial bundle
# and an audited executable. No solver or initial model is loaded from disk.
from salmonpy.tddft import TDDFTInitialState
with TDDFTInitialState.from_salmon(
        'existing-initial-bundle', 'rt.inp', executable=engine,
        max_steps=1000, threads=1) as model:
    with model.new_state(dt=0.02) as resumed:
        resumed.restore_checkpoint('checkpoint-0005')
        current = resumed.step([0.1, 0., 0.])
```

`save_checkpoint(directory)` and `restore_checkpoint(directory, max_bytes=None)`
are TDDFTState-only methods returning None. Saving refuses existing destinations
and paths inside the original input bundle or either owned initial bundle. Disk restore allows a separately
constructed model with identical contents. Exact package version, input, all
initial-file hashes, executable hash, dt, initial A, capacity and thread count
must match; rejection preserves the current worker and accepted state.

Storage is Zarr v3 with uncompressed little-endian float64 arrays, time-chunked
A/current histories and a final manifest.json hashing metadata and every chunk.
No Zarr Python dependency is required. Initial GS/potential payloads are not
copied: only initial filenames and hashes are recorded. Retain the existing
initial bundle separately. Stored current is internal matter current J;
the public current property continues to return charge current j.

Ordinary save failures clean up the new destination. Incomplete writes left by
forced termination are rejected. Power-loss durability, concurrent writers and
external edits during loading are not promised. The default byte ceiling is
derived from the receiving model capacity: 48 bytes per step plus edge-chunk
padding, fixed vectors, metadata and hash inventory. It is not a machine-wide
fixed number. Optional max_bytes also applies the caller budget to stored files
and decoded array payload. Required files/directories follow the format; extra
entries are rejected. Neither bound is a process RSS limit: additional copies, two replay
workers and GS memory are excluded. Restore still replays the entire accepted history with
existing current checks, so its cost grows with saved step count. Replay failure
preserves the original worker.

Same-machine single-thread synthetic regressions check byte-identical 20-step
continuation from another Python process; these are interface tests, not
converged-GS physical validation. Coupled hosts must also rewind field states.

## Explicit exclusions

This is not the handoff's Zarr GS archive, an output NPZ, a cross-version
restart, a C ABI or a variable-dt method. Do not
serialize/edit tokens or pass output archives to `restore`. Snapshot public
metadata (`kind`, `schema`, `version`, `step_index`, `time`) is diagnostic only;
private fields are not an editable checkpoint schema.

Embedded TDDFT, native E/H FDTD and MPI sessions do not gain this API. A host manually coupling supported material/field states must save and
restore **both at the same accepted coupling boundary**. There is no automatic
transaction across components; a single state's atomic restore does not provide
MPI/global rollback. `couple()` still creates and closes its own states.
Changing dt, adding adaptive stepping/subcycles or making native TDDFT reentrant
requires separate work. Transverse prescribed-A use remains the basic workflow.

SBE now keeps native model/density handles; snapshot export and atomic restore
retain the same contract. See [SBE ownership](SBE.en.md) for close/GC and memory tradeoffs.

LD retains native pole histories and transfers them only for snapshots/restoration. See [LD.en.md](LD.en.md) for ownership and memory.
