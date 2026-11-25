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

## Near field from prescribed dipole sources

Select the SALMON input dipole sources with `from_salmon(..., native_near_field=True)`.
The default is disabled. The source settings are `num_dipole_source`,
`cood_dipole_source`, `vec_dipole_source` and `rad_dipole_source`.
A validation record for this executable must explicitly include near-field comparisons;
ordinary uniform-field evidence does not enable the feature.
Short C2H2 tests verified batch dipoles, independent states, live restoration and direct disk restart after process exit. These are software tests, not physical convergence checks.

The scope is PZ, fixed ions, scalar pseudopotentials and one rank.
At most two sources are supported, with a positive finite radius and `Acos2` or `Ecos2` pulses.
Source coordinates and radius follow the SALMON input unit system; Python `step(E)` remains in atomic units.
The input spatial profile and pulse envelopes are retained. `step(E)` supplies an **additional uniform field**.
Zero sources use the uniform-field path. Restoration requires matching source inputs and mode.

This exposes prescribed dipole sources, separately from general spatial E/H input.

## Spatial E/H near-field mode

`from_salmon(..., near_field=True)` couples spatial E/H to molecular electrons.
The current scope is **PZ, fixed ions, unpolarized scalar pseudopotentials, CPU, one rank**.
Magnetic coupling is orbital. Spin/SOC magnetic coupling and radiation feedback to FDTD are outside this interface.
The molecular region is vacuum: Gaussian atomic units give `B=H`.

### Direct FDTD grid input without an expansion

`GridFieldSeries` directly interpolates saved E/H and evaluates the electric path integral of
Takenaka, Taketsugu and Iwasa. No multipole degree or decay rate is selected. The current provider
uses a uniform3D Cartesian grid with trilinear interpolation, splitting integration at component-grid
cell faces. Orbital magnetic coupling uses the same interpolated fields.

```python
from salmonpy.near_field import GridFieldSeries

# Collocated SampledField objects. H_left is at n-1/2, H_right at n+1/2.
grid = GridFieldSeries.from_yee_samples(
    E_samples, H_left, H_right,
    counts=(16, 17, 17), origin=(7., -15., -15.), spacing=(2., 2., 2.),
    molecule_position=(21.2, .7, .9), center=(0., 0., 0.), time_origin=.303)
grid.write('NEW_INPUT_BUNDLE/salmonpy-near-eh.dat')  # preparation only, reuse a saved GS

with IsolatedTDDFTInitialState.from_salmon(
        'NEW_INPUT_BUNDLE', 'rt.inp', runtime=runtime,
        max_steps=20, near_field=True) as initial:
    with initial.new_state(dt=.003) as state:
        for _ in range(7):
            state.step([0., 0., 0.])  # additional uniform E; grid fields come from prepared history
        state.save_direct_checkpoint('NEW_CHECKPOINT')
```

Coordinates/spacing are bohr, time is atomic time. Counts are sample-node counts, not the entire
FDTD simulation size. Origin is the first sample position. molecule_position is the FDTD position
of the integration center; center is that same point in TDDFT coordinates. Fortran maps
`FDTD position = molecule_position + (TDDFT position-center)`; neither position must be a grid node. Axes must share their orientation; rotations are not inferred.
The patch must cover the electronic box, stencil halo and nonlocal-projector paths, not just the atoms.
No extrapolation or implicit periodic wrapping occurs. The molecular region is vacuum with H=B.

If H is already centered at E times, use `from_samples(E_samples,H_samples,...)`.
Structured arrays use `GridFieldSeries(origin,spacing,times,E,H,...)`, with E/H shaped
`[time,z,y,x,component]`. Staggered spatial origins use `[E/H,component,axis]`; extract a common
patch with the same counts/spacing for all components. Structured half-step H arrays can use
`from_yee_arrays(...,H_left_times=...,H_right_times=...)`. Temporal centering, interpolation and
integration run in Fortran; missing H or endpoints are not inferred.

Input size scales with grid points times frames; extract only the needed local patch and time window.
Provide relative times starting at0 with constant dt matching TDDFT and frames0..max_steps+1. Update n reads
frames n and n+1 and advances through `step(E_extra)`. Do not use the coefficient-only `step_near_field`.
The original field time is retained as time_origin in the input file; it does not reproduce prior molecular
exposure. Specify the initial electronic and illumination conditions. Direct checkpoints retain grid data,
coordinate mapping and electronic state; restore with the same model, dt and calculation program.

Use `python -m salmonpy.examples.prepare_grid_field --help` for preparation only; no solver is launched.
The executable requires direct-grid validation evidence. Orbital magnetic coupling is tested separately
from the original electric-only paper; spin coupling and feedback radiation are outside this endpoint.
Interpolation cannot reconstruct unresolved FDTD information; check field-grid and molecular convergence.

Analytic staggered affine fields/path integrals agree within about2.1e-17 absolute error.
For a smooth field, spacing0.5/0.25/0.125bohr gives relative RMS E errors6.60/1.46/0.328% and
potential errors3.42/0.646/0.0901%. Actual saved LD/FDTD fields driving an off-grid C2H2 molecule
match continuous execution after a process-exit restart at step7 through20 for all nine observables
and output streams. These verify connection/restart, not physical convergence of a real nanostructure.

### Choose the expansion center during fitting

```python
from salmonpy.near_field import NearFieldSeries

# Matching-position/time SampledField objects: E:au/H:au, positions in bohr.
near = NearFieldSeries.from_samples(
    E_samples, H_samples, center=(0., 0., 0.), radius=8., lmax=1)
print(near.center, near.relative_residual)  # columns: E, H
near.write('GS_BUNDLE/salmonpy-near-eh.dat')  # preparation only; no TDDFT launch
```

An explicit molecular center is recommended. If omitted, Fortran chooses the sample centroid in the first fit;
if radius is omitted, it uses the farthest sample distance. Center and radius remain fixed for all frames and restarts.
Radius normalizes the basis; it does not improve the underlying FDTD resolution.

Each Cartesian E/H component includes real solid harmonics with radial terms.
`lmax` is the maximum total degree d. Terms are ordered d=0..lmax, then angular
l=d,d-2,..., then m=-l..l, multiplying `(r/radius)^d Y_lm`. Thus d=l+2k includes the
radial factors r^(2k). There are P=(lmax+1)(lmax+2)(lmax+3)/6 terms per component;
`E_coefficients`/`H_coefficients` have shape `[time,P,xyz]`.
Their units differ from the electric-potential-only coefficients. lmax=0 is uniform; lmax=1 includes field gradients,
thus electric quadrupole coupling. Uniform H produces orbital magnetic-dipole coupling.
Higher orders are available, but a finite expansion need not represent an arbitrary field exactly.
Residuals above the default 5% warn while still returning usable fits.

For Yee half-step H, call `NearFieldSeries.from_yee_samples(E_samples,H_left,H_right,...)`.
E is at n, H_left at n-1/2, H_right at n+1/2. The positions/times are checked and H is averaged in Fortran.
Supply both initial brackets explicitly; missing endpoints are not guessed or extrapolated.


Molecular positions and expansion centers are real coordinates; they need not coincide with FDTD grid points.
Fit surrounding E/H samples about the chosen center without moving the molecule to a grid point or
substituting a nearest-neighbor field. This is a one-way conversion of saved fields into TDDFT input.
The electric path-current map and checked FDTD injection are described below; the orbital magnetic return interface is described below.
Samples should cover the molecular electronic-density region. A single smooth expansion across a metal
interface, or extrapolation outside the sampled region, is not an accuracy guarantee.

The required order for tips and bow ties has not been established. Vary apex radius, gap, molecular
position and fitting region as well as order; check absolute and relative errors at held-out points.
The radial-complete finite basis spans all polynomials through the requested total degree.
It does not recover information unresolved by FDTD or enforce Maxwell consistency. Check residuals
at independent positions and the FDTD resolution before claiming molecular accuracy.


The earlier diagnostic compared the angular-only representation with complete polynomials.
The TDDFT input now includes radial terms through the same total degree. Small grid-point
residuals still do not establish between-grid or molecular-response accuracy. Insufficient sample
count or rank is rejected as unavailable, not treated as convergence. At degree8, P=165 independent
coefficients per component require at least165 full-rank sample points. The V3 field file and V7
direct checkpoint identify this representation; older field files/checkpoints must be regenerated.

A connection test reused the C2H2 GS for20 updates with nonzero quadratic radial E/H terms.
After saving update7 and terminating the process, all nine observables for the remaining13 updates
matched continuous execution exactly, with zero additional evolution during restoration.
Independent analytic fields of total degree2/4/6/8 gave off-sample E/H errors at most about1.3e-15 au.
These verify representation and connection, not physical convergence for tips or bow ties.


### Fields with directional exponential decay

```python
near = NearFieldSeries.from_samples(
    E_samples, H_samples, center=(0., 0., 0.), radius=8., lmax=4,
    decay_vector=(0., 0., .25))  # inverse bohr; decay along positive z
near.write('GS_BUNDLE/salmonpy-near-eh.dat')  # preparation only
```

Each E/H component uses `exp(-decay_vector dot (r-center)) (|r-center|/radius)^d Y_lm`.
Radius normalizes the polynomial part; the decay vector stays in inverse bohr. The default
zero vector retains the polynomial representation. E and H share the vector, fixed across frames
and restarts along with center/radius/degree. Pass the same vector to `NearFieldFrame.from_samples`.
No automatic rate selection occurs. If a rate is selected using separate tuning points, retain
independent final validation points; the selected rate need not identify a unique physical decay length.
One factor is not a general tip/bow-tie field model, and grows in the opposite direction; evaluate
within the validated region covering the molecule.

Fortran separates field representation from molecular coupling integrals. Electric coupling integrates
E along the center-to-electron path; magnetic coupling integrates the path-parameter-weighted H.
Stable exponential moments avoid small-argument cancellation; magnetic projector phases use an
exponential-aware line quadrature. A common sampled-path integration entry is available for future
representations; this does not implement direct FDTD interpolation or two-way feedback.


Custom field representations can call `salmonpy.near_field.couple_path_samples(displacement, parameters, weights, E, H)`.
Displacement is in bohr; parameters run from0 to1; weights integrate that dimensionless parameter.
E/H are collocated `[points,3]` arrays on the path. Results are electron scalar potential energy and A/c.
The caller supplies field samples and quadrature nodes/weights and checks quadrature convergence.
This function does not advance TDDFT.

Preparation example: `python -m salmonpy.examples.fit_near_field --help`.
The execution example `examples/isolated_response.py` accepts `--near-field`.
Fitting, writing input and importing never launch TDDFT.

### Pass coefficients at each update

```python
with IsolatedTDDFTInitialState.from_salmon(
        'GS_BUNDLE', 'rt.inp', runtime=runtime, max_steps=80,
        near_field=True) as initial:
    with initial.new_state(dt=.003) as state:
        saved = state.savepoint()
        frame = state.step_near_field(near.frame(1), near.frame(2))
        state.restore_savepoint(saved)
        state.release_savepoint(saved)
        state.step_near_field(near.frame(1), near.frame(2))
        state.save_direct_checkpoint('NEW_CHECKPOINT')
```

Update n consumes field frames n and n+1 in SALMON's existing order.
Frames carry center, radius, order, time, step index and E/H coefficients.
`NearFieldFrame.from_samples(...,step_index=n,center=...,radius=...)` fits one spatial E/H frame;
your own code may also construct a `NearFieldFrame` directly from coefficients.
To revise a shared field frame already used by a preceding update, restore a savepoint from before that update.

For an unknown future field, use the initial frame's `write_live_input(path,max_steps=...,dt=...)`
and call `step_near_field` for every update. This held initial template is not a forecast.
Ordinary `step(E_extra)` uses the already supplied E/H history plus the additional uniform E argument.

Use live savepoints or direct disk checkpoints. E-only replay snapshots do not encode live spatial coefficients
and are not offered in near-field mode. The coefficient history, center, radius, magnetic arrays and nonlocal
phases are saved directly. Restoration does not repeat electronic evolution.

### Coupling and validation scope

For displacement r from the chosen center, Fortran evaluates
`V_e=r·integral_0^1 E(s r) ds` and `A/c=-r×integral_0^1 s H(s r) ds/c`.
V_e is the electron scalar energy. The symmetrized kinetic operator retains orbital magnetic coupling and A²;
nonlocal projectors use line-integral phases from their ion centers.
This is a prescribed classical-field [radial-gauge potential construction](https://arxiv.org/abs/physics/0204034).

Short C2H2 tests cover electric quadrupole/magnetic effects, live versus prescribed coefficients,
restoration to the initial state and direct restart after process exit. They do not establish material-response
convergence. Check expansion order, FDTD/electronic grids, dt and fit residuals for the target observable. Independently fitting E/H does not automatically enforce Maxwell consistency.
A field-free GS exposed to nonzero initial E/H evolves from that specified initial state;
preceding irradiation history or a field-dressed GS is not created automatically.

## Quasistatic electric-potential approximation of saved fields

`salmonpy.near_field.ElectricPotentialSeries` fits collocated E samples using real regular solid harmonics in Fortran.
**It approximates the scalar electric-potential component, rather than providing general near-field E/H coupling.**
It does not include H coupling, magnetic multipoles, radiation feedback or information unresolved by the FDTD grid.

```python
from salmonpy.near_field import ElectricPotentialSeries

# electric_samples: SampledField('E', times, positions, values, 'E:au')
# Positions are in bohr; align component locations and times first.
potential = ElectricPotentialSeries.from_samples(
    electric_samples, center=(0., 0., 0.), radius=8., lmax=2)
print(potential.relative_residual)
potential.write('GS_BUNDLE/salmonpy-ylm.dat')  # preparation only, no solver launch
```

An explicit `reference=vacuum_samples` subtracts matching-position/time vacuum samples with the same incident source;
subtraction also runs in Fortran. Relative residuals above the default 5% warn and remain available.
Near a zero field, also compare `electric_field(...)` with the original samples.
Increasing order does not generally restore missing magnetic or rotational field components.

The real basis uses the [NIST complex normalization](https://dlmf.nist.gov/14.30#E1):
positive m uses `sqrt(2)*(-1)^m*Re(Ylm)`; negative m uses the corresponding positive-order imaginary part with the same factor.
Order is l=1..lmax, then m=-l..l; lmax=1..8 is supported.
`V_e=sum(c_lm*(r/radius)^l*Y_lm)`, `E=grad(V_e)`; coefficients are Hartree.
The constant l=0 term is omitted. The regular expansion assumes no external charges inside the region.

Execution is a separate explicit operation and requires spherical-harmonic validation evidence for this executable.

```python
with IsolatedTDDFTInitialState.from_salmon(
        'GS_BUNDLE', 'rt.inp', runtime=runtime, max_steps=80,
        electric_potential=True) as initial:
    with initial.new_state(dt=.003) as state:
        for _ in range(40):
            state.step([0., 0., 0.])  # additional uniform E
        state.save_direct_checkpoint('NEW_CHECKPOINT')
```

Coefficient times start at zero with constant dt and cover indices0..max_steps+1.
Update n reads coefficients n and n+1 in SALMON's existing sequence.
Use `restore_direct_checkpoint` with matching model, dt and executable.
For a selected FDTD time window, retain the original field step/time and the molecular time origin explicitly.
Starting a field-free GS at that window does not reproduce a molecule already exposed to its preceding history.

A 27-site 3D LD/FDTD saved-field test reached80 molecular updates. Direct restart after process exit at update40 matched
continuous execution exactly. Uniform-field harmonics matched the batch reference within about4e-12 e*bohr.
However, the actual FDTD tests had large field-approximation residuals; these are not near-field physical-accuracy validations.
Use the spatial E/H near-field mode above when both electric and magnetic coupling are needed.

## Distributing charge response for spatial electric fields

`GridFieldSeries.charge_polarization` and `charge_current` implement the transpose of the electric path-integral operator in Fortran. They return polarization or interval-average current on the E-component grids. Input `charge_density` is **signed physical charge density**, e/bohr³, not positive electron number density. Electrons alone have negative charge density. Supply point quadrature volumes in bohr³ through `volume_elements`. Fixed ionic charges cancel from the current difference only for unchanged ionic positions, sample points and volumes.

```python
J = field_series.charge_current(
    molecular_positions, old_charge_density, new_charge_density,
    volume_elements=volumes, dt=dt)
P = field_series.charge_polarization(
    molecular_positions, new_charge_density, volume_elements=volumes)
```

Paths translate TDDFT coordinates relative to the expansion center to the FDTD molecule position and split at component-grid cell faces. Output shape is `[z,y,x,component]`; P has units e/bohr², J atomic units. Components use the input E grids. The spatial work identity is `dt*sum_grid(E*J*dV)=sum_molecule(delta_rho_charge*V_e*dv)`, where `V_e(r)=(r-R).integral(E(R+s*(r-R)),s=0..1)` is the electric path integral entering the electron Hamiltonian. This identity for supplied arrays does not establish temporal accuracy or full-domain energy balance.

This is an array-level calculation component, usable with charge densities supplied by custom code. Explicit worker-density acquisition and direct component-grid current injection are now available below. Live grid input and electric-only own-vacuum feedback are validated in the development driver below; a general production coordinator is not provided. Do not send this J through DipoleMap or the macro-site current entry: that would apply another interpolation. The separate orbital magnetization-current interface is described below; complete coupled charge-continuity validation remains open. Use `python -m salmonpy.examples.grid_electric_response --help` for saved-array preparation.

## Read TDDFT charge density and pass component currents

`state.charge_density()` explicitly reads the current accepted electronic charge density, without propagation, reinitialization or replay. Its owned, read-only `ChargeDensitySamples` contains `step_index`, `time`, `positions` in bohr, negative electronic `charge_density` in e/bohr³ and `volume_elements` in bohr³. Ions are excluded; fixed ions cancel from endpoint differences.

```python
old = state.charge_density()
state.step(E)
new = state.charge_density()
J = series.charge_current(new.positions, old.charge_density, new.charge_density,
                          volume_elements=new.volume_elements, dt=state.dt)
```

This is supported for validated molecular PZ, fixed-ion, single-rank executables. The public loader verifies the executable SHA256 and `charge_density` evidence. It is not automatic per-step output: each explicit request updates a binary transfer file inside the private worker directory. Account for this I/O at scale. Permanent density storage remains opt-in. Direct electronic checkpoints already contain density; readout after restoration uses that accepted step without evolution.

`field.advance_edges(indices,current)` directly consumes physical current density already placed on E components. `indices[N,3]` are zero-based global grid offsets; `current[N,3]` is in atomic units, without another volume or sign factor. Each field rank consumes only owned cells, so list lengths may differ. Duplicate rows add. This entry currently requires checked mode and retains the existing recovery gates and E->LD->H update. Do not add a second macro-site interpolation or pass a current defined on a different component grid.

`electric_edge_geometry()` returns global counts, origins[component,axis] in bohr and spacing. `electric_edges(indices)` explicitly reads E components and sums owned values within the field group only; no application-wide MPI communication is added. Current is averaged over the interval about to be advanced. A coupling scheme may require trial/rollback and an own-vacuum reference: these interfaces alone are not a self-consistent spatial feedback driver.

Real C2H2(18,432 samples) readout left observations bit-exact against an independent control; process-exit restart at electronic step2 matched densities and observations exactly. Charge changes were converted in Fortran and radiated through one/two field ranks. Direct injection matched the existing deposition; NaN rejection, direct live rollback/retry and a separate-process direct field checkpoint continuation matched array controls. These are software connection tests, not physical convergence or complete spatial molecular feedback validation.

## Streamed grid input and explicit step numbers

A runtime with matching validation evidence accepts `state.step_grid_field(window,step_index=n)`. `window` is a two-frame GridFieldSeries with relative times `[0,dt]` and `time_origin=n*dt`. Its E/H frames represent FDTD steps **n,n+1**; the call advances molecular step **n-1 to n**. Origins, spacing, counts, molecule position and integration center must match the prepared input.

```python
window = GridFieldSeries(origin, spacing, [0., dt],
                         [E_n, E_next], [H_n, H_next],
                         molecule_position=position, center=center,
                         time_origin=n*dt)
state.step_grid_field(window, step_index=n)
```

Fortran validates the target step, fixed geometry/dt and every value before installing either frame. Wrong geometry, NaN/Infinity and incomplete transfers preserve numerical state. Electric-only input must cover electronic integration paths. Before first enabling H, stencil-halo and nonlocal projector paths must also fit; insufficient coverage rejects H while electric-only continuation remains available. Orbital magnetic input is supported; this does not add spin/SOC magnetic response.

Initialization still requires salmonpy-near-eh.dat. For unknown future fields, explicitly prepare max_steps+2 zero or otherwise chosen placeholder frames with the fixed geometry, then replace each required pair. Allocated history remains proportional to max_steps*patch_cells. Use a molecular patch and estimate memory/private transfer I/O for long runs.

The existing isolated_response example accepts `--grid-series FIELDS.npz`. Supply origin, spacing, times, E, H, molecule_position and center with frames0..N+1. E_NPZ retains N additional uniform-field vectors and dt. Use --direct-checkpoint / --resume-direct; E-only propagation replay cannot reconstruct spatial input. Step labels are local to this run. Recording a source-time origin does not recreate earlier illumination of the GS.

## Validated electric spatial feedback

`tools/test_grid_feedback_mpi.py` combines public density readout, streamed input and direct current injection in a finite-scope validation driver. Total and own-source-only vacuum FDTD groups are separate. Fortran subtracts their fields and maps molecular charge changes to polarization current; only owning field ranks consume their current rows. Full grids are not broadcast to every application rank.

First hold FDTD once with zero current. Molecular update0->1 consumes FDTD frames1,2 and finishes at field2. Update95->96 consumes frames96,97 and finishes at field97. A trial field predicts the drive; restoring fields and injecting the response must yield a matching selected drive. Current support is placed in vacuum. These tests pass zero H to the molecule: this is **electric-only feedback**, not complete magnetic-return E/H coupling.

Off-grid C2H2 passed outer MPI4/6 checks. The192-molecular-step LD case produced an environment-return field difference about4.0e-7 au; predictor/accepted selected-field disagreement was about1.4e-17 au. Rejection/retry and process-exit restoration at molecular96/field97 matched continuous arrays and ordinary output streams bit-exactly, with0extra evolution on restoration. This validates software connections, not physical convergence or arbitrary material layouts. Unexpected errors terminate this validation MPI job; it is not a general production recovery coordinator.

## Return spinless orbital magnetic response

The electric interaction remains in the length gauge: a uniform field gives the electron scalar potential `r·E`. With magnetic fields, the magnetic vector potential remains in the kinetic term `(p + A_mag/c)^2/2`, where `A_mag = -r × ∫₀¹ s B(R+s r) ds`. This does not switch the electric drive back to the velocity gauge. A uniform magnetic field gives `A_mag = B × r/2`; conversely `B = curl A_mag`. Recovering the supplied field through that identity requires a solenoidal field and allowance for interpolation error.

With matching executable-bound evidence, `state.orbital_magnetic_response()` reads the occupied-state orbital magnetization without propagation. Initial scope: spinless fixed-ion molecule, orthogonal cell, direct E/H grids and one electronic worker rank. Spin magnetization is excluded.

Fortran computes `M = -δ〈H〉/δB` from the actual applied Hamiltonian: local paramagnetic and diamagnetic terms and the derivative of the nonlocal projector path phases. Input H denotes vacuum B at the molecule. In Gaussian/atomic units the returned magnetic charge current is `j_M = c curl M`.

```python
previous_magnetic = state.orbital_magnetic_response()  # molecular step n-1
state.step_grid_field(window, step_index=n)
magnetic = state.orbital_magnetic_response()
# electric_current contains ONLY the density-derived electric polarization current
current = magnetic.total_current(electric_current, previous=previous_magnetic)
# Route these E-component samples to their owning field ranks, then advance_edges once.
```

`magnetization[z,y,x,3]` uses the supplied H-component grids. `current_density()` returns the instantaneous magnetic current alone. With `previous=...`, adjacent accepted magnetizations are averaged in Fortran to obtain the interval current. `total_current(...,previous=...)` adds that interval current once to the electric polarization current. Do not supply an already complete orbital current: that would double count the magnetic part.

Current conversion requires Yee origins: E shifted half a cell along its own component and H along the other two axes. Include padding around magnetization support and its circulating current; magnetization must vanish at positive patch faces. A truncated curl that would discard surface current is rejected. Off-grid molecule positions are supported. Spatial derivatives, path adjoints, interval averaging and current summation run in Fortran.

Magnetization is derived from saved wavefunctions and applied fields; no extra physical history is introduced. After direct restoration, read the previous response immediately, then combine with the next accepted response. C2H2 molecular steps2->4 reproduce magnetization, interval currents and ordinary observables bit-exactly without restoration evolution.

`tools/test_orbital_magnetic_response.py` checks finite differences of the actual frozen-state Hamiltonian and direct restart. `tools/test_electric_edges_mpi.py response --response-history ...` injects the resulting electric/magnetic/combined currents into FDTD. These are connection checks. The default `test_grid_feedback_mpi.py` mode is electric-only; its `--orbital-magnetic` own-vacuum extension is described below. General failure recovery and physical convergence remain separate work.

## Bidirectional coupling with orbital magnetic response and an own-vacuum reference

`tools/test_grid_feedback_mpi.py --orbital-magnetic` connects direct E/H input and orbital magnetic-current return to both total and own-vacuum FDTD groups. Verified scope: spinless fixed-ion PZ C2H2, one electronic worker and one/two ranks in each field group. Both groups receive the same molecular current. Subtracting their E/H removes the vacuum self-field and retains the environment-scattered return. This policy does not include vacuum radiation reaction.

`field.component_field_geometry()` returns distinct E/H component origins. `field.component_fields(indices)` returns `[2,N,3]`: raw E and H averaged from adjacent half-step brackets at the same integer field step. Spatial staggering is preserved. Fortran retains the preceding raw H independently of observation buffers; rollback and disk checkpoints include it. A single owning-field-group reduction collects E/H; full grids are not replicated across all application ranks.

First hold both FDTD groups once with zero current. Molecular update n-1->n uses field frames n,n+1. A zero-current trial predicts frame n+1; Fortran obtains electric polarization current and the interval orbital magnetic current, and adds them exactly once. Restore both fields, inject the same combined current into each group's owned E components and accept only when the selected E/H matches the prediction.

The192-molecular-step trial saves molecular96/field97, including both FDTD/LD states and the molecule. A fresh MPI launch restores them and reads the previous magnetization from the restored quantum state, without extra evolution. E/H, molecular observables and six ordinary streams matched continuous execution bit-exactly. Rejection/retry at step64 and one/two field ranks per group also matched. Omit `--orbital-magnetic` for the existing electric-only case.

```sh
mpiexec -n 6 python tools/test_grid_feedback_mpi.py feedback \
  FIELD_BUILD E_WORKER VALIDATION_JSON EXISTING_GS NEW_OUTPUT \
  --orbital-magnetic --ld --field-ranks 2 --steps 192 --save-step 96
```

Resume with matching inputs, a new output directory and `--resume SAVED_OUTPUT/checkpoint`. The executable-bound receipt needs `grid_eh`, `grid_stream`, `charge_density` and `orbital_magnetic` evidence. Installed wheels place the tool in the package's `_runtime_source/tools/`. This uses a fixed validation geometry/source; it is not a general arbitrary-geometry or failure-atomic production coordinator. Unexpected MPI errors terminate the test job. Physical convergence, spin magnetization and vacuum radiation reaction remain separate work.

## One-way spatial FDTD E/H excitation

Use `python -m salmonpy.examples.molecular_grid_response WORKER VALIDATION BUNDLE input.inp FIELDS.npz NEW_OUTPUT`. This reads a saved FDTD history and evolves only the molecule. It does not start FDTD, return molecular current, subtract a self-vacuum reference or calculate radiation. Molecular back-action on the illuminating structure is neglected.

Prepare the existing GS bundle with matching `salmonpy-near-eh.dat` grid geometry and a worker validation record declaring `grid_eh` and `grid_stream`. FIELDS.npz contains `origin`, `spacing`, `times`, `E`, `H`, `molecule_position`, `center` and optional `time_origin`. E/H shape is `[frames,z,y,x,3]`; distinct Yee-component origins are preserved. Supply time-centered H from the FDTD component-field reader, not uncentered half-step H. Atomic units and vacuum at the molecule apply.

Frames are numbered 0..N+1 at constant dt. Molecular step n-1 -> n uses frames n,n+1. `history.window(step_index=n,dt=dt)` constructs that pair without fitting, spatial averaging or temporal interpolation. The Fortran worker performs spatial interpolation and the molecular path-integral coupling. Source `time_origin` is recorded as provenance; it does not recreate earlier molecular excitation.

Use `--stop-after 2 --direct-checkpoint NEW_CHECKPOINT`, then run the same command with a new output directory and `--resume-direct NEW_CHECKPOINT`. Restore requires the same full history, prepared geometry, dt, capacity and runtime. The saved state retains its staged fields; resumption starts at the next molecular step without replay. Inspect `restore_statistics` and `final_statistics` in the output metadata. Outputs include molecular step numbers, dipole, electronic change, electron number and, when available at every frame, internal energy. Availability is recorded. No auxiliary uniform-E archive is required.
