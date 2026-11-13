# Individual research tasks with SALMONpy

All new numerical inputs are atomic units. These Python helpers are included in
wheels; only standalone LD requires the bundled native library. Full FDTD uses
an explicitly supplied, separately built SALMON executable.

## Geometry and FDTD

```python
from salmonpy.geometry import Grid, Shape, Geometry
from salmonpy.fields import AcosPulse
from salmonpy.fdtd import FDTD, Dielectric, PlaneWave, SalmonExecutable

grid = Grid(counts=(32, 32, 32), spacing=(2., 2., 2.))
geometry = Geometry(grid, (Shape('ellipsoid', size=(8., 8., 8.)),))
material_ids = geometry.material_ids()   # independent material-grid preview
xy = geometry.slice(axis=2, index=16)   # usable with a plotting library
pulse = AcosPulse(amplitude=.001, omega=10., duration=.05)
calculation = FDTD(
    geometry, materials=(Dielectric(epsilon=2.),), dt=.001, steps=200,
    source=PlaneWave(pulse, axis=0, position=-12.),
)
print(calculation.input_text())          # review/export without running
# engine = SalmonExecutable('/absolute/path/to/salmon', validated_sha256)
# result = engine.run(calculation, '/new/output/directory')
# table = result.observations[0]
# Ez = table.column('E_z')
```

`Shape.size` is full diameter/width, not radius. Coordinates are isolated SALMON
material-grid centers; staggered E/H component locations differ. Later shapes
replace earlier material IDs. Only axis-aligned boxes/ellipsoids are exposed.
The preview was compared to SALMON's emitted cube files, including mixed odd/even
noncubic grids. Use material IDs 1..len(materials); 0 is vacuum.

`Boundaries` specifies six faces in x-,x+,y-,y+,z-,z+ order, each PML or PEC.
The legacy solver's PML is fixed at eight cells; clearance includes the interface
plane that receives its coefficients. Material/PML overlap is rejected because
this initial interface assumes vacuum PML. Periodic structures, rotated geometry,
nonvacuum PML, oblique incidence and adjustable PML thickness are not yet exposed.

FDTD accepts an `Impulse` with LD media, or a `PlaneWave` Acos2 pulse with LD or
constant `Dielectric` media. The actual SALMON impulse mode does not support
constant dielectric objects, so this is rejected before execution. Plane sources
must be transverse, inside the nonabsorbing region, and clear of matter. Positive
x/y/z propagation axes only. Source planes are snapped by the legacy solver to
its grid; inspect stdout.log for resolved coordinates.

The common pulse follows TDDFT's Acos2 convention. FDTD's existing cos-carrier
convention differs by a quarter cycle; the adapter adds 0.25 cycles to its CEP so
the same Python pulse represents the same continuous waveform. Kernels are not
changed. `amplitude` is an envelope scale, not an exact peak electric field.

Model, source, grid, dt and observations are fixed in the immutable FDTD object.
Execution-only choices (new output directory, timeout and thread count) belong to
`run`. Every run records settings, binary/input SHA256 and status. Existing output
directories are never overwritten. Hashes identify a binary, not establish its
safety: source/object audits belong to its build/validation. This is a batch
interface, not an in-process FDTD step API or an MPI launcher.

## Pulse creation independent of a solver

```python
import numpy as np
from salmonpy.fields import AcosPulse, sample_pulses
pulses = [AcosPulse(.01, .2, 100., start=10.),
          AcosPulse(.01, .2, 100., start=130.)]
A = sample_pulses(pulses, np.arange(0., 250., .1))
E = A.electric_field()
# A.write_salmon('new-pulse.data')  # four columns for ae_shape='read', unit_system='au'
```

A means SALMON `Ac=A/c`, so E=-dAc/dt. The discrete electric field uses backward
differences at A.time[1:], as SALMON does; it is an interval-average field assigned
to the new endpoint. No initial E is invented. CEP is in cycles. Pulse arrays
are read-only copies. SALMON's file interpolation returns zero outside the file's
range (including its final endpoint); make exported A zero at the endpoints.
Composed pulse files are not yet accepted by the FDTD wrapper.

## Output and reproducible spectra

`analysis.read_table(path)` handles numbered SALMON column headers, retained
unit labels and Fortran D exponents. It makes no silent unit conversions.

```python
from salmonpy.analysis import read_table, spectrum
# table = read_table('obs1_at_point_rt.data')
# spec = spectrum(table.column('Time'), table.column('E_z'), window='hann',
#                 time_unit=table.units[0], value_unit=table.units[3])
```

The transform is dt*sum(w*x*exp(+i omega t)), including the actual time origin.
Only nonnegative frequencies are returned. There is no automatic DC subtraction,
window gain correction, one-sided doubling, or conversion to HHG/emission power.
`spec.config` records the window, sign, normalization, grid and units. Irregular
sampling is rejected rather than silently resampled.

## A-only crystalline backend contract

`crystal.AInitialState` and `AResponseState` are typing protocols for replaceable
backends. `initial.new_state(dt=..., time0=..., initial_A=...)` must create an
independent state without modifying the reusable initial condition. The state
owns wavefunctions/density matrices, history, previous A and current, and exposes
`step(A_new)`, `current`, `time`. Only A is an input; E is not an alternative.
Backend coordinates, gauge and J volume normalization must be converted to the
common laboratory-frame atomic-unit convention at the adapter boundary.

`crystal.respond(initial, A_samples)` creates a new state and calls that same step
method for each sample. It validates uniform times, current shape/finiteness and
backend time advancement. Accumulated floating-point time roundoff is tolerated
with a bound capped at 0.001 time steps; a materially wrong step is rejected. It preserves the backend's initial current rather than
assuming it is zero. A backend must guarantee independent states and atomic step
failures; Python cannot enforce those promises for arbitrary external engines.

SALMON TDDFT/SBE adapters and initial-file loaders are **not yet implemented**.
The tests use explicit interface test doubles, not physical TDDFT/SBE models.
The next extraction must retain native initial-state formats, MPI ownership,
wavefunction/density-matrix histories, rotations, gauge and time alignment. DFT
initial-state generation remains a separate stage; never advertise a fake solver
or silently substitute LD for electronic dynamics.

## Built-in velocity-gauge SBE

```python
import numpy as np
from salmonpy.sbe import SBEInitialState
from salmonpy.crystal import respond
from salmonpy.fields import AcosPulse

initial = SBEInitialState.from_salmon(
    'gs_export', 'Si', bands=32, kpoints=512, electrons=32,
    lattice=np.diag([10.26]*3), active_bands=32,
    correction_order=0, nonlocal_correction=False)
drive = AcosPulse(amplitude=0.01, omega=0.056, duration=100.).sample(
    np.arange(2001)*0.05)
result = respond(initial, drive)
# Independent repeat: no excitation history is carried over.
repeat = respond(initial, drive)
a = initial.new_state(dt=0.05, initial_A=drive.values[0])
b = initial.new_state(dt=0.05, initial_A=drive.values[0])
J = a.step(drive.values[1])  # b remains at its initial state
```

Lattice vectors are columns, all quantities are in atomic units, and `A` is
SALMON's `Ac=A/c`. `density` and `current` return immutable snapshots.
The loader reads single-spin, orthogonal-format SALMON text exports with
`esp[a.u.]`; eV and nonorthogonal-format files are rejected. Array construction
accepts any right-handed lattice and uses `(band,k)` energies and
`(band,band,xyz,k)` matrix elements. Occupations are exactly SALMON SBE's closed
shell rule: two electrons in each of the first `electrons/2` bands, irrespective
of the export's occupation column. `active_bands` defaults to all file bands.

The shared Fortran kernel preserves four sequential Taylor terms at midpoint A,
then computes endpoint J, including correction orders 0–3 and optional nonlocal
velocity corrections. This is not an exactly unitary propagator; no adaptive
step or norm correction is introduced. `dt` must be checked for convergence.
The initial output is the physical current at initial A, not a forced zero.

Python SBE is serial, with independent per-state arrays and no global native
calculation state. Each step currently copies initial arrays into a temporary
Fortran context; this first adapter prioritizes equivalence over large-system
throughput. A single evolving state should not be stepped concurrently.
The native TDDFT adapter, length gauge and restart/checkpoint loading remain
unimplemented. This loader creates a fresh ground-state density, not a restart.
