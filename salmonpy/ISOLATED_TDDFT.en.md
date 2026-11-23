# Isolated TDDFT: E input and named dipoles

Explicit process-isolated PZ, fixed ions, unpolarized scalar pseudopotentials and middlepoint propagation. Crystalline TDDFT/SBE remain A-driven. The response is a charge dipole, not a current density. Grid, dt and electronic-box convergence remain the user's responsibility.

## Calculation program and initial state

```python
from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState

runtime = load_isolated_runtime('salmon', 'validation/report.json')
with IsolatedTDDFTInitialState.from_salmon(
        'existing-bundle', 'rt.inp', runtime=runtime,
        max_steps=100, threads=1, timeout=60) as initial:
    with initial.new_state(dt=.02) as state:
        saved = state.savepoint()
        frame = state.step([0., 0., 1.e-5])
        print(frame.time, frame.dipole, frame.electron_number)
        state.restore_savepoint(saved)  # same worker, no evolution calls
        state.release_savepoint(saved)
```

E/time use atomic units and dipoles use e*bohr. `dipole` includes ions and electrons; `electronic_change` is the electronic dipole change. Arrays are read-only. Initial E/time are zero. `new_state(position=...)` declares the field evaluation position (default zero, bohr), without moving or rotating atoms.

Loading checks a bounded local validation record and the executable SHA256. Current schema: `salmonpy.isolated-worker-validation/5`. E3 still returns nine numerical values, but evidence of same-process direct restoration with zero extra evolution calls and stale-token rejection is required. Revision3 receipts are rejected. Loading never downloads, builds, runs validation, initializes MPI or starts a worker. A receipt is local evidence, not a signature or a substitute for eigensolver auditing.

`tools/test_isolated_worker.py EXECUTABLE EXISTING_REFERENCE NEW_VALIDATION` explicitly runs ordinary response and direct-restore validation, reusing the supplied GS without SCF. It is a development calculation tool. This does not announce a distributed runtime catalog: supply a locally validated calculation program.

The bundle needs an isolated RT input, `data_for_restart/info.bin`, `occupation.bin`, `wfn.bin` and relative pseudopotential assets. The loader owns an immutable initial copy; states own independent working copies/processes. Existing states survive closing their loader. No periodic-GS conversion or implicit isolated SCF is performed.

## Direct savepoint in the same worker

`savepoint()` retains both orbital buffers, density/potential histories, Poisson work state, dipole/field histories, energy and accepted index in Fortran and returns `IsolatedSavepoint`. One live token per state. `restore_savepoint(token)` restores the numerical state and output journal in the same process, without propagation. It keeps the token for repeated restoration. Explicitly call `release_savepoint(token)` when done; release alone preserves the current state.

Foreign, released or pre-restart tokens are rejected before restoration. Tokens are not a disk format. Release a live token before using replay `snapshot()/restore()` or `save_checkpoint()`. The caller owns these lifetimes; closing the state ends its savepoint too.

`worker_statistics()` explicitly queries PID, real evolution count `updates`, accepted index `step` and `direct_restores`. No application-wide communication is added to normal updates. Restore alone does not increment `updates`.

While a direct token is active, failures never start a replay worker. A temporary output-journal failure suspends updates (`usable=False`); explicitly finish restoration with the same token. Worker loss or an indeterminate restore is reported as unusable. Do not resume updates without checking usability. This is not arbitrary OS/MPI failure recovery.

Only the two audited principal output streams are currently supported. Unsupported spatial output modes reject savepoint creation only; ordinary E-driven calculations remain usable. The backup needs memory for an additional electronic state.

## Common field input and observations

```python
from salmonpy.field_input import FieldPacket, FieldValue
frame = state.step_fields(FieldPacket(E=FieldValue(
    [0., 0., 1.e-5], time=state.time+state.dt,
    position=(0.,0.,0.), unit='E:au')))
```

Run this inside a live state context. E time/position/units are checked; missing E differs from zero. A/H are not consumed and E is not implicitly reconstructed from A. Standalone `step(E)` retains the existing constant-interval E integrator order.

`frame.electron_number` is the actual Fortran density integral. Relative loss1e-5 is an advisory diagnostic, not an automatic stop. `frame.internal_energy` is present only at fresh endpoints according to the existing energy cadence, otherwise `None`. Never silently add energy evaluations or subtract zero-field drift.

## Direct disk restart after process exit

```python
with model.new_state(dt=dt) as state:
    for E in first_fields:
        state.step(E)
    state.save_direct_checkpoint('new-electronic-checkpoint')
# The original worker has now exited. Reopen the same GS model/settings.
with model.new_state(dt=dt) as resumed:
    resumed.restore_direct_checkpoint('new-electronic-checkpoint')
    assert resumed.worker_statistics()['updates'] == 0
    for E in remaining_fields:
        frame = resumed.step(E)
```

`salmonpy.isolated-direct/1` stores the Fortran numerical buffers/histories,
Poisson and XC work state, accepted index, cached observations and output streams.
Per-file SHA256, model settings and executable identity are checked before loading.
The original GS bundle is supplied separately; it is loaded for initialization,
with no SCF and no replay of elapsed electronic steps. A fresh candidate process
is installed only after complete import and output reconnection. Failed imports
leave the current worker unchanged. Imported states never silently recover by
replay. One additional electronic state may coexist during restoration.

Use the same executable, GS/input, dt, capacity, position and thread settings.
This binary format has no cross-compiler or changed-layout portability claim.
Release live savepoints before disk operations. `max_bytes` optionally limits
bytes read, not resident memory. Saved directories are never overwritten.
Only the existing audited output modes are supported. The current worker receipt
revision5 includes real process-exit restart proof; older receipts are rejected.

The response example accepts `--stop-after 7 --direct-checkpoint NEW_CHECKPOINT`.
A second invocation with the same E NPZ and settings can use
`--resume-direct NEW_CHECKPOINT`; it advances only the remaining samples.
These operations save the molecule only; a complete coupled-field restart must
also coordinate both FDTD states and the coupling clock.

## History persistence and replay (separate from direct restore)

`snapshot()/restore()` save and replay accepted E/observation history. They require the same initial model, dt, position, capacity, threads and binary. A candidate is checked before replacement, preserving the original if replay fails. Cost grows with history, with two workers coexisting during validation. Ordinary `step()` recovery without a live direct token also uses this history-replay mechanism.

```python
state.save_checkpoint('new-history-checkpoint')
resumed.restore_checkpoint('new-history-checkpoint', max_bytes=16*1024*1024)
```

The caller creates `resumed` with the same initial data/settings. `salmonpy.isolated-tddft-checkpoint/3` is a hashed Zarr v3-compatible history format containing E0, initial/current nine observations, E history and observation history. Supply the GS separately. `max_bytes` adds a file-reading budget, not an RSS limit. This is not a direct wavefunction checkpoint. Use the direct API above for numerical persistence after worker exit.

## Outputs and FDTD

`python -m salmonpy.examples.isolated_response EXECUTABLE VALIDATION_JSON BUNDLE RT_INPUT E_NPZ NEW_OUTPUT` writes a dipole archive from E[steps,3] and scalar dt in atomic units. `--checkpoint` explicitly saves replay history. Retain units/provenance and store results outside the source tree.

Use [Molecular/FDTD coupling](DIPOLE_RADIATION.en.md), `MolecularFeedback` and the complete `examples/molecular_feedback.py` for feedback. Cancellation uses `savepoint/restore_savepoint/release_savepoint` without electron-history recomputation. Own-vacuum subtraction also removes vacuum radiation reaction.

The short C2H2 PZ MPI6 test uses64 field intervals: normal execution needs63 electronic updates; one rejected interval and retry needs64. The PID stays constant within each run. Restore adds zero evolution calls and every accepted saved array/six-face record exactly matches the previous result. This is wiring validation, not converged optical accuracy.
