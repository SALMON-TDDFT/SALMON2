# Isolated TDDFT: E input and named dipoles

`salmonpy.isolated_tddft` provides an explicit process-isolated PZ, fixed-ion,
unpolarized scalar-pseudopotential, midpoint interface. Periodic crystals remain
A-driven through `salmonpy.tddft`. This interface returns **charge dipoles**, not
current density. Grid/time/box convergence is the user's responsibility.

## Runtime and initial state

```python
from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState
runtime = load_isolated_runtime('/path/to/salmon', '/path/to/validation/report.json')
with IsolatedTDDFTInitialState.from_salmon(
        '/path/to/bundle', 'rt.inp', runtime=runtime, max_steps=100,
        threads=1, timeout=60) as initial:
    with initial.new_state(dt=.02) as state:
        saved = state.snapshot()
        frame = state.step([0., 0., 1.e-5])
        print(frame.time, frame.dipole, frame.electronic_change, frame.electron_number)
        state.restore(saved)
```

All values are atomic units: E in atomic electric-field units, time in atomic
time, dipoles in e*bohr. `dipole` is total charge dipole (ions plus electrons);
`electronic_change` follows the solver electronic dipole-change convention.
Readonly arrays and named fields prevent confusing these with periodic matter
current. Initial E and time are zero. Spatial position defaults to (0,0,0) bohr
and can be set with new_state(position=...); it labels the uniform field sample,
not a translation/rotation of molecular coordinates.

Runtime loading checks an explicit local validation record bounded to 64 KiB, its E3
protocol/capabilities, successful comparison/isolation/replay evidence, and exact
binary SHA256. It performs no download, compilation, MPI initialization or worker
launch. validation records are local evidence, not cryptographic certification or a source
code audit. Obtain them by running `tools/test_isolated_worker.py` from a source
checkout against an explicitly built/audited standard runtime and the small
reference fixture. validation records are specific to that binary and tested scope, not to
the accuracy of arbitrary inputs. The schema is
salmonpy.isolated-worker-validation/1; ordinary A-only runtime manifests are not
sufficient. No prebuilt E-enabled catalog asset is published by this change.

The bundle contains an isolated RT input, data_for_restart/info.bin,
occupation.bin, wfn.bin and its relative pseudopotential/other assets. Use solver
records compatible with the executable. Loading copies a private immutable
baseline; each independent calculation states it and owns a writable working copy and child process.
External edits do not affect the captured baseline. The solver worker rejects
unsupported input modes during initialization before propagation. The loader
does not convert periodic GS records or generate an isolated GS automatically.

## Common field packet and restoration

```python
from salmonpy.field_input import FieldPacket, FieldValue
# Inside the live state context above.
frame = state.step_fields(FieldPacket(E=FieldValue(
    [0., 0., 1.e-5], time=state.time+state.dt,
    position=(0.,0.,0.), unit='E:au')))
```

E must match the declared time, position and units. Missing is not zero. A/H are
not consumed by this calculation implementation. `step(E)` is the explicit canonical-unit shorthand.
The E value is used at the solver midpoint propagator's step index, preserving
its update order; it is not inferred by integrating a supplied A.

`snapshot()` returns an in-memory token of accepted E/dipole history. `restore()`
requires the same initial-model object, dt, position, capacity, threads and binary.
It creates a candidate process and replays history; only after all observations
match (rtol=1e-12, atol=1e-14) does it replace the old state. Failed restoration
leaves the original state untouched. Candidate and old worker coexist during
restore, so budget for two workers. States can restore after the initial loader
is closed, while its leased baseline remains alive. Closing the loader prevents
new states but does not stop existing ones.

A failed step attempts to reconstruct the previous accepted state and raises
RuntimeFailure; if recovery fails the state closes. Input validation does not
advance it. FieldAdapter itself still invalidates on calculation implementation exceptions, so
rebuild an outer adapter if needed. There is no solver wavefunction snapshot, MPI whole-session restoration, or multi-context Fortran
process guarantee here. Replay cost increases with history length.

## Archive example and FDTD

```sh
python -m salmonpy.examples.isolated_response EXECUTABLE VALIDATION_JSON BUNDLE RT_INPUT E_NPZ NEW_OUTPUT
```

E_NPZ contains real E[steps,3] at dt,2dt,... and scalar positive dt, in atomic
units. The archive includes time, E (including zero at t0), named dipoles and units.
Use frame.dipole with DipoleMap for one-way radiation; see DIPOLE_RADIATION.en.md.
Do not apply the periodic current sign reversal or divide by the TDDFT vacuum
volume. Self-field feedback is a separate model choice, not enabled automatically.

## Disk replay checkpoint

```python
# On a live state; use a new directory outside the initial bundle.
state.save_checkpoint('new-isolated-checkpoint')
# A new state from a separately reconstructed identical initial model:
resumed.restore_checkpoint('new-isolated-checkpoint', max_bytes=16*1024*1024)
```

`salmonpy.isolated-tddft-checkpoint/3` stores E0, observables0, observables, E_history
and observable_history in the strict uncompressed Zarr v3 subset. Observation vectors
have nine components: electronic_change3, total charge dipole3 (e*bohr), measured
electron number1, internal energy1 (Hartree), and freshness1.
It is a replay record, not a wavefunction checkpoint or postprocessing archive.
GS files and the validated runtime must be supplied separately. Initial file
hashes, executable hash, dt, capacity, threads, position and package version must
match; initial-model object identity may differ. A periodic A checkpoint is rejected.
Hash/layout/size checks precede candidate replay; failure leaves the live state
unchanged. max_bytes bounds archive/decoded payload, not process RSS. Default
limits follow configured capacity. Only new destinations are accepted.

The response sample accepts `--checkpoint NEW_DIRECTORY` to save its final state.
To continue later, recreate the same capacity/settings, restore, then supply the
remaining E steps. A sample saved at its capacity can only replay up to that end;
plan extra capacity in your own run when further propagation is required.


### Measured electron number

E3 exposes `frame.electron_number` and `state.electron_number`: the integral of
occupied electron density, not an estimate from input nelec. Reuse the result
already calculated by Fortran dipole observation; no additional global reduction.
Replies contain nine values: electronic dipole change3, total charge dipole3,
electron number1, internal energy1 and energy freshness1. Arrays are read-only and the scalar must be finite.
Evidence uses `salmonpy.isolated-worker-validation/3`; checkpoints use
`salmonpy.isolated-tddft-checkpoint/3`. Observation arrays have six dipole columns
followed by measured electron number, internal energy and a 0/1 freshness flag. Old E1/checkpoints are rejected.
Electron number participates in replay/restore comparisons. The molecular FDTD
fixture diagnoses loss in Fortran with target1e-5; exceeding the target alone does
not stop a calculation. Nonfinite values reject the trial.


### Internal energy and work

`frame.internal_energy` and `state.internal_energy` return internal energy in
Hartree, excluding external-field interaction energy, only when the existing
endpoint eigen-energy evaluation is fresh. Other endpoints return `None`.
Never silently change out_rt_energy_step or force additional evaluations.

The development molecular FDTD fixture accumulates the actual held-endpoint work
Σ E[n]·(p[n]−p[n−1]) and midpoint-field work separately in Fortran. Compare internal
energy only at fresh endpoints. A zero-field drift control is measured separately
and is not silently subtracted. Finite residuals are diagnostics, not stop gates.
In NPZ output, internal_energy_raw needs energy_fresh; alternatively select
energy_observation_steps/energy_observations. This is distinct from full Maxwell
energy balance. Old E1/E2 and old checkpoint formats are rejected.

The molecular FDTD fixture defaults to midpoint-E coupling. Its Fortran average supplies the constant interval E to unchanged standalone step(E). Endpoint driving is an explicit numerical control. See DIPOLE_RADIATION.en.md for clocks, work and own-field scope.
