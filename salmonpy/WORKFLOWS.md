# Individual research tasks with SALMONpy

All new numerical inputs are atomic units. These Python helpers are included in
wheels; standalone LD and built-in SBE require the bundled native library. Full FDTD uses
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
# engine = SalmonExecutable(executable_path, validated_sha256)
# result = engine.run(calculation, '/new/output/directory')
# table = result.observations[0]
# Ez = table.column('E_z')
```

`Shape.size` is full diameter/width, not radius. Coordinates are isolated SALMON
material-grid centers; staggered E/H component locations differ. Later shapes
replace earlier material IDs. Only axis-aligned boxes/ellipsoids are exposed.
The preview was compared to SALMON's emitted cube files, including mixed odd/even
noncubic grids. Use material IDs 1..len(materials); 0 is vacuum.

`Boundaries` specifies six faces in x-,x+,y-,y+,z-,z+ order: PML, PEC, or paired periodic faces. See [transverse-periodic slabs](FDTD_MATERIALS.en.md).
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
interface, not an in-process FDTD step API. Existing native MPI FDTD can be launched with FDTDMPI; see FDTD_MATERIALS.en.md.

## Pulse creation independent of a solver

```python
import numpy as np
from salmonpy.fields import AcosPulse, sample_pulses
pulses = [AcosPulse(.01, .2, 100., start=10.),
          AcosPulse(.01, .2, 100., start=130.)]
A = sample_pulses(pulses, np.arange(0., 250., .1))
E = A.electric_field()
# A.write_salmon('new-pulse.data')  # four columns for ae_shape1='input', unit_system='au'
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

The built-in velocity-gauge SBE adapter and text loader are described below.
The process-isolated TDDFT adapter is described below. Common-interface tests also use an
explicit test double; SBE has separate native reference tests.
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

Python SBE uses one cached native read-only model and independent persistent
state/trial density buffers. Normal steps exchange small vectors; density and
snapshots explicitly export arrays. Use state.close or a context manager; GC
also releases abandoned states. See SBE.md/SBE.en.md for ownership and memory
tradeoffs. A single evolving state should not be stepped concurrently.
SBE length gauge and excited-state restart/checkpoint loading remain
unimplemented. This loader creates a fresh ground-state density, not a restart.

## Native TDDFT with isolated state

`TDDFTInitialState` implements the same `new_state`/`step(Anew)` contract as SBE.
It keeps one persistent SALMON process per evolving state because the complete
solver still owns module globals, communication scratch space and XC histories.
The existing `initialization_rt` and `time_evolution_step` execute unchanged inside
that process. This is not a TDDFT shared library bundled in the wheel.

Prepare a self-contained directory with a valid periodic RT input, pseudopotential
files and `data_for_restart/{info.bin,occupation.bin,wfn.bin,...}` from GS. Use
relative paths for all assets and a plain `sysname`. The initial-state loader
copies this bundle once and records hashes. States lease that shared baseline;
each worker receives an independent writable copy whose complete manifest is
verified before launch. Changes to the original files do not change the loaded initial state. Binary
records are read by SALMON itself and retain its compiler/format constraints.
GS and RT records share formats: the caller must supply a GS bundle; the loader
does not infer ground-state convergence or automatically classify arbitrary
wavefunctions as GS versus excited states.

```python
from salmonpy.fdtd import SalmonExecutable
from salmonpy.tddft import TDDFTInitialState
from salmonpy.crystal import respond

# Use the SHA256 of your source/object-audited, worker-enabled executable.
executable = SalmonExecutable(executable_path, validated_sha256)
with TDDFTInitialState.from_salmon(
    'prepared_bundle', 'rt.inp', executable=executable,
    max_steps=len(drive.time)-1, threads=1, timeout=300.) as initial:
    result = respond(initial, drive)  # closes its process, also on failure
    with initial.new_state(dt=0.05) as a, initial.new_state(dt=0.05) as b:
        Ja = a.step([0.001, 0., 0.])
        Jb = b.step([0.002, 0., 0.])  # no state or native globals shared with a
```

Current supported scope: periodic, unpolarized, fixed ions, **PZ, TB-mBJ, Libxc-PZ or Libxc-PBE**, CPU, one MPI
rank with optional OpenMP, `middlepoint`, transverse supplied total A/c, fresh GS
rather than an excited RT checkpoint. The existing predictor-corrector option is
retained. `libxc_pz` and `libxc_pbe` require a separately built SALMON with
`USE_LIBXC=ON`; use the same functional and pseudopotential for GS and RT.
The loader does not certify GS provenance or detect a GS/RT functional mismatch.
Other functionals, AETRS/HSE, longitudinal/internal-Maxwell feedback,
spin-orbit, +U, projection/unfolding and moving ions are rejected. `time0` must
be zero; `dt`, capacity and initial A are fixed when creating a state. Capacity
controls native history-array allocation and cannot be enlarged mid-run.

The input's pulse, base directory and restart output settings are replaced only
in worker mode; the public drive is exclusively supplied endpoint A. The existing
midpoint Hamiltonian, predictor-corrector, density/potential/history updates,
and endpoint-current calculation retain their order, including existing
nonlocal-phase conventions. Since 2026-11-14 the state returns the electric
(**charge**) current density j, as does SBE; SALMON's electronic matter current
J = −j stays internal (`state.matter_current`). No numerical integrator or
unitarity correction is introduced.

Always close direct states, preferably with `with`. `respond` does this for you.
Invalid A leaves the state unchanged. On worker failure the adapter reconstructs
and replays accepted A samples before raising; if recovery fails it closes the
state explicitly. Replay is expensive and is not a scientific checkpoint format.
Snapshots and per-process copies also require disk space. N live states retain
N+1 managed bundles (one shared baseline plus N worker copies), excluding the
caller input and retained output archives. Closing the initial loader refuses new
states but existing leases keep replay available. The baseline is removed after
the loader and final state release it. During restore an extra worker copy exists;
the replay candidate shares the same baseline. No symlinks/hardlinks are used.

For native TB-mBJ use `xc='tbmbj'`; Libxc is not required. The native
`cval` policy is preserved (positive fixed value, otherwise density-derived).
Keep this policy fixed in the initial bundle. Native TB-mBJ uses predictor-corrector
even if `yn_predictor_corrector='n'`. On the tested coarse Si grid, default
Gaussian initial orbitals produced NaNs; a converged PZ wavefunction used as an
initial guess followed by full TB-mBJ SCF convergence succeeded. The worker must
load the resulting TB-mBJ GS, not merely relabel a PZ GS. See
`docs/tbmbj-validation-2026-11-13.md`.

## Postprocessing outputs

See OUTPUTS.md for `respond(..., output_directory=...)`, `state.finish(directory)`,
selected native spatial outputs and the versioned `OutputArchive` readers.
Native TDDFT save failure retains a recovery directory rather than discarding
completed outputs. Existing run directories are never overwritten.

バンド計算・描画は [BANDS.md](BANDS.md) を参照してください。

Optional, independent ASE structure/path/band-result bridge: [ASE.md](ASE.md).

Prepare a converged GS from ASE or plain arrays: [DFT.md](DFT.md).

Optional ASE Calculator with native PZ/Libxc-PZ/Libxc-PBE forces (KY/NC-UPF) and opt-in wavefunction warm start: [ASE_CALCULATOR.md](ASE_CALCULATOR.md).

Field components with the j→A protocol (SALMON tr / lo / 2d and 1D multiscale propagation) and the `couple` driver: [FIELDS.en.md](FIELDS.en.md). Interfaces exchange the charge current j; components declare their time layout (salmonpy.time/1), checked by `couple`.

Experimental host-MPI TDDFT and native E/H coupling now have a separate entry: [COUPLED.en.md](COUPLED.en.md). One point per solver group and one native lifetime per process; do not infer unrestricted reentrancy or restart support.

## ホスト向けAPIと状態保存

[公開API・例外](API.md)と[メモリ上のsnapshot／restore](SNAPSHOTS.md)を参照してください。
0.1.0a5からLD・SBE・4種類のPythonの場で固定dtの保存・復元を使えます。
ネイティブFDTD・TDDFT・MPIの再開を追加するものではなく、電流・時刻規約も変更しません。

Persistent standalone LD ownership, fixed-dt coefficients and close/context-manager usage: [LD.en.md](LD.en.md).

Portable caller-supplied GS arrays (GS/2), separate from native GS/1; explicit paired Gamma converter supports serial TDDFT initialization: [GS_ARCHIVE.en.md](GS_ARCHIVE.en.md).

DFT-generated serial Gamma GS/1 export: GroundStateResult.export_archive with a new paired exporter; real/complex orbitals retained, density reconstructed and provenance recorded. See GS_ARCHIVE.en.md.

SerialRestartConverter opts into aligned orthogonal multi-k native/portable GS conversion with explicit verified k-point files and weighted density; GammaRestartConverter stays Gamma-only. See GS_ARCHIVE.en.md for regular-source vs explicit-import scope and current tolerances.

### Serial TB-mBJ ground states

Use `DFT(..., xc="tbmbj", cval=None)` or `dft_from_atoms(..., xc="tbmbj", cval=1.2)`. Material JSON recipes accept xc and optional cval. None selects the automatic coefficient; a finite value greater than 1e-10 selects a fixed coefficient. The coefficient is retained in GS records, warm-start compatibility checks and TDDFT inputs. TB-mBJ is a model potential: total energies and forces are unvalidated diagnostics, with a warning rather than a run prohibition. ASE Calculator optimization and native parallel SCF support are unchanged. Run `python -m salmonpy.examples.dft_tbmbj --help` for the serial example. Validation here covers synthetic native XC evaluation and interface contracts, not real-material accuracy.

## Optical conductivity and dielectric response: one polarization

Use `analysis.optical_conductivity(time, charge_current_x, total_E_x,
current_convention="charge", field_kind="total", units="atomic_gaussian",
window="hann", drive_rtol=1e-8, metadata={"material_id":"Si","polarization":"x"})`.
Then call `analysis.dielectric_function(sigma)`. Results have omega, values,
valid and config. Inputs are simultaneous scalar time series; A is not E.
For A-driven crystals the caller must provide E with consistent differentiation
and time placement. Convert SALMON matter current J explicitly using j=-J.
Only Gaussian atomic units are accepted. One selected drive gives a single
polarization response, not a complete tensor; one circularly polarized run does
not determine all tensor elements.

The existing compiled NumPy FFT uses exp(+i omega t), normalized by dt. With
exp(-i omega t), sigma=j/E and epsilon=1+4*pi*i*sigma/omega. Fortran computes
ratios and diagnostics. Drive bins at or below drive_rtol times the maximum,
and nonfinite conversion results, carry NaN with valid=False. Dielectric DC is
undefined. An all-zero drive yields masked output rather than stopping a solver.
External-field ratios are not silently treated as bulk dielectric functions.
The same window is applied to both signals; window, time range, dt, peak E and
caller provenance are retained. Windowed finite-time ratios are estimates;
check duration/window dependence, especially for broadband drives.

`optical_linearity(reference, candidate, rtol=...)` compares supplied normalized
conductivities on common valid frequencies, reporting relative RMS. The caller
must ensure matching material, polarization and drive shape. Zero reference
response has undefined relative error/acceptance. `optical_sum_rule(sigma,
expected_weight=..., normalization="...")` integrates Re sigma between adjacent
valid bins with the trapezoidal rule. For an all-electron Gaussian-AU system,
the positive-frequency f-sum reference is pi*n/2. Pseudopotentials, nonlocal
current operators and missing bands require an appropriate caller-specified
reference. Masked gaps and unobserved frequency ranges are not extrapolated.
Neither diagnostic is a physical accuracy certificate or launches another run.

Run `python -m salmonpy.examples.optical_response` for synthetic data, not a
converged real-material optical calculation.

LD accepts optional epsilon_infinity (default1), retained by material JSON/providers and checkpoint identity. P/J remain dispersive poles; displacement and background_energy explicitly include the instantaneous background. See LD.en.md for field coupling and optical-current conventions.
