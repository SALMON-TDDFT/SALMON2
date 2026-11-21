# Custom FDTD materials and interfaces

[日本語](FDTD_MATERIALS.md) · Topic: `fdtd_materials` · Paired revision: `2026-11-21.1`

Run `python -m salmonpy.examples.fdtd_layers MANAGED_RUNTIME NEW_DIRECTORY`
with a verified runtime v2. This one-process, one-thread example contains two
synthetic finite layers. It is not a converged infinite-film optical calculation.

Define each isotropic dispersive material with `LorentzDrude(plasma_frequency=...,
strength=[...], resonance=[...], damping=[...])`. All quantities use atomic units;
frequency and damping are angular frequencies. Each material may have its own
number of poles. A zero resonance describes a Drude pole. Use validated material
coefficients; explicit fitting is available in LD/OPTICAL_DATA, without a convergence guarantee.
For a constant isotropic medium use `Dielectric(epsilon=..., mu=..., conductivity=...)`.
See [LD conventions](README.en.md) before translating coefficients from another
code. The example coefficients are deliberately synthetic.

Pass the materials as a tuple to `FDTD`. Shape ID 1 refers to the first material,
ID 2 to the second; ID 0 is vacuum. Boxes and ellipsoids can be combined. If IDs
are renumbered, change both the tuple and shapes consistently. At overlapping
shape cells the last shape wins: this geometrical order is intentional. Material
IDs do not encode physical priority. Inspect the generated material-ID cube.
Keep material clear of the vacuum PML, and test grid/time convergence. Current
limits are 200 materials, 100 poles per LD material, and 1000 shapes.

The default interface rule is a deterministic staircase assignment, without averaging
unlike dispersive models:

- An electric-field edge between two nonvacuum cells takes the lower-coordinate
  cell's complete material model. If either cell is vacuum, retain the legacy
  vacuum-edge coefficients.
- A magnetic-field site takes the first nonvacuum cell in fixed spatial order.
  For Hx this is `(i,j,k)`, `(i,j+1,k)`, `(i,j,k+1)`, `(i,j+1,k+1)`; Hy and Hz
  use the corresponding transverse planes. A fully vacuum stencil stays vacuum.
- Native MPI initialization exchanges faces, edges and corners of the material
  grid. The existing E/H and LD time-update order remains unchanged.

This rule preserves the tested single-material/vacuum calculations. Results at
interfaces between unlike materials intentionally change from the old numeric-ID
priority. It does not assert subpixel accuracy, rotational invariance or improved
convergence order. At observation points the legacy printed material parameters
refer to adjacent **cells**, not effective interpolated E-edge coefficients.

Runtime v1 is rejected for multi-material runs. Obtain/build v2 explicitly; imports
never compile. `salmonpy-run.json` records `interface_assignment`. A new `run`
creates independent solver state. Recreate the calculation when changing material
parameters or timestep.

This public batch interface supports LD and constant dielectric materials. It does
not yet accept arbitrary current callbacks, anisotropic tensors, tabulated optical
response, periodic boundaries, oblique incidence or TDDFT/SBE materials in E/H
FDTD. The Python batch launcher can also start the existing native MPI FDTD with
`FDTDMPI`; see below. TDDFT/SBE A-to-J interfaces and the native A-form multiscale path remain
separate capabilities.


## MPI FDTD with many local LD states

FDTD was already MPI-parallel in SALMON. This API selects its existing spatial
partition; it does not add a new Maxwell parallel algorithm. Each FDTD rank owns
its local field grid and LD P/J/Q histories, selecting coefficients by material
ID. Several LD materials and many LD cells can reside on one rank. No dedicated
MPI group is allocated per LD cell. Native implicit E/LD/H update order is kept;
the separate endpoint-A LD response-service adapter is not used here.

```python
from salmonpy.fdtd import FDTDMPI, SalmonExecutable
solver = SalmonExecutable.from_mpi_build("MPI_BUILD_DIRECTORY/build.json")
result = solver.run(calculation, "NEW_OUTPUT", mpi=FDTDMPI((2, 2, 1)), threads=1)
```

Run this from a fresh serial Python controller, not under mpiexec or initialized
mpi4py. FDTDMPI's three integers define the spatial process grid; their product is
the MPI rank count. `launcher=("mpiexec",)` can be replaced with an argument tuple
for the host launcher; there is no shell interpolation. Each rank needs at least
two cells per direction. MPI launching currently targets macOS/Linux.

Build explicitly using `python -m salmonpy.runtime build` with USE_MPI=ON and the
host MPI C/Fortran compilers. See RUNTIME.en.md; choose the build job count from
your memory/CPU budget. The successful build record binds MPI capability to the
binary hash. An older record without that capability requires a fresh build;
editing CMakeCache does not confer MPI capability. Portable serial runtime bundles
are not MPI builds. For externally audited builds outside this builder, the
low-level SalmonExecutable(..., mpi_enabled=True) is the caller's explicit
capability assertion, not automatic binary detection.

`salmonpy-run.json` records the launcher arguments, rank grid/count, threads and
completion/failure/cancellation. Timeout and interruption terminate the local
launcher process group; remote job cleanup still depends on the MPI/scheduler
launcher. Do not share one output directory between jobs. Threads is the OpenMP
setting; account for any separately configured BLAS threads in CPU planning.

```sh
python -m salmonpy.examples.fdtd_mpi_layers MPI_BUILD_JSON NEW_DIRECTORY --grid 2 2 1
```

This example uses three synthetic LD materials with 1, 2 and 3 poles. It prints
material-cell counts and leaves memory explicitly unmeasured. Production sizes
need resource planning and numerical convergence tests.

For electronic feedback into E/H, use the separate [native coupling entry](COUPLED.en.md). Do not put TDDFT objects in the batch API material list.

## Transverse-periodic slabs

`Boundaries(('pml','pml','periodic','periodic','periodic','periodic'))` uses x PML and y/z periodic boundaries. Periodic faces must come in opposite pairs and use zero-phase wrapping. Coordinates and material-ID previews stay centered. Fortran field communication selects each axis independently without changing electronic periodicity. Existing vacuum-PML material checks remain. Use a runtime with this per-axis implementation rather than an older bundle.

`python -m salmonpy.examples.fdtd_periodic_slab RUNTIME NEW_OUTPUT` selects a dielectric; add `--ld` for LD. Explicit `tools/test_fdtd_periodic.py AUDITED_MPI_BUILD_JSON NEW_RESULTS` compares vacuum-normalized complex Fresnel transmission. The solver remains Fortran; Fourier sums and the analytic reference are evaluated in test-only C. On the 256×8×8 fixture with 1 and 2 ranks, maximum relative complex error over 3 frequencies was about 0.0554% for the dielectric and 0.0375% for LD. Transverse probes and MPI 1/2 results agree. This verifies the stated synthetic slab, not arbitrary materials/conditions. The common grid initializer still has a minimum-domain constraint on thin MPI decompositions; this example uses 8 transverse cells.

Set an LD instantaneous background with `LorentzDrude(..., epsilon_infinity=2.25)`; default1 preserves prior results. Non-unit-background Impulse propagation continues, while legacy normalized spectra are omitted with a warning. See [LD.en.md](LD.en.md) for PlaneWave propagation, whole-domain CFL and pole P/J versus D conventions.

Fitted LD valid_frequency_range is retained. Out-of-band source carriers and broadband impulses warn without blocking propagation. Supply an E spectrum to LD.check_excitation for a fuller spectral-band check.

## Common isotropic tensor configuration

`model.tensor_config` / `config.tensor_config` returns LDTensorModelConfig with one independent
component, epsilon_ij=epsilon*delta_ij. `LorentzDrude.from_config(tensor)`, FDTD materials
and LDMaterialProvider accept it. This applies across modes, not only3D E/H.
`permittivity_tensor(omega)` returns immutable complex(...,3,3) arrays; scalar APIs remain.
FDTD.tensor_materials and its run record expose the representation. Native FDTD now supports isotropic, diagonal and reciprocal off-diagonal tensors; standalone/MPI tensor material response is now implemented below.
See OPTICAL_DATA.en.md for oxide data, source assumptions and limits.

## Diagonal tensor FDTD and waveplates

FDTD accepts TensorLorentzDrude whose principal axes align with grid xyz.
Fortran derives each E/LD-current semi-implicit coefficient from that axis's
background and oscillator weights. P/J/Q update ordering and Yee interface
ownership are unchanged. Reciprocal off-diagonal propagation is also connected below; no component is silently zeroed.

The explicit input is media_type='ld-diagonal' plus salmonpy-ld-diagonal.dat.
When exporting input_text(), export input_assets() too. SalmonExecutable.run and
NativeEHField stage assets automatically; their identity is included by SHA256.
Obtain/build a matching runtime explicitly; old runtimes are rejected before use.
Scalar epsilon_em/f_ld contain first-axis compatibility values; the asset holds
all independent diagonal values. Isotropic input has no extra asset and retains
the prior input/arithmetic path.

```sh
# Preview only; passing this does not launch a calculation on a management node.
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --plate quarter
# Explicitly start inside an allocated compute job
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --plate quarter --run
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON OTHER_DIRECTORY --plate half --run
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON REFINED_DIRECTORY --spacing 4 --run
```

The sample is a synthetic lossless slab with n_x=1.5,n_y=2, not real-material
optical-accuracy validation. analysis.waveplate_diagnostic uses Fortran to extract
carrier phasors and normalized Stokes from an E trajectory containing the full
pulse/ringdown. Its reference is normal incidence on a vacuum-bounded plane slab
with equal in-phase transverse incident amplitudes, including Fresnel and multiple
reflections, not just geometrical retardance. S3=2 Im(Ex conj(Ey))/S0; phasors use
exp(+i omega t). One-frequency agreement is not a convergence certificate.
Nonunit/diagonal LD impulse calculations retain propagation and raw data while
omitting the legacy scalar optical normalization.

The reported phasor is a trapezoidal Fourier integral with units E(au)*time(au), not the pulse peak amplitude. Polarization uses its normalized ratios.

## Reciprocal off-diagonal native E/H FDTD

Pass TensorLorentzDrude with lab-frame symmetric SPD epsilon_infinity and PSD
pole strengths. Each pole keeps a scalar resonance/damping shared by its matrix
components. This is not antisymmetric magneto-optics. It applies to the native E/H
backend with 1D-like, 2D-like or 3D geometry; it does not add an anisotropic background
to the separate maxwell.py A-only field components. Native nonzero Bloch phase is
not implemented for the new tensor operator (Python FDTD does not expose it).

Full tensors use media_type='ld-tensor' and salmonpy-ld-tensor.dat, staged by both
batch and NativeEHField. Matching runtime capability and asset SHA256 are checked.
All six independent entries are retained. On staggered Yee edges, the material
operator uses paired edge-to-node means and transpose means, projected onto each
material's own edges. A diagonal-preconditioned conjugate-gradient solve couples
E; the same pole operator drives J before the existing P and Q updates. Numerical
work is Fortran. Diagonal/isotropic input retains the existing arithmetic path.

FDTD(..., tensor_rtol=1e-10, tensor_maxiter=200) controls this solve. A finite solve
that exhausts its budget warns and continues; inspect its residual and refine as
needed. Nonfinite arithmetic or nonpositive curvature rejects the trial, with
native host state rollback permitting retry. PEC fixed edges stay excluded.
The full-tensor direct checkpoint binds the immutable tensor model and solver
settings; mismatches are rejected before committing. Old scalar/diagonal checkpoint
streams retain their format. This is not repartitionable field checkpointing.

```sh
# Preview only; no solver or output directory is started
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --angle-degrees 30
# Explicit run in allocated resources; incident principal-axis amplitudes are equal
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --angle-degrees 30 --run
```

waveplate_diagnostic(..., angle=pi/6) accepts lab E; its phasor and Stokes results
are in lab xy and its reported phase difference is in the principal axes. angle is
in radians; the sample CLI uses degrees. Complete the pulse/ringdown time window.
This is a synthetic Fresnel-slab validation, not real-material certification.

Local regression (not required CI):
```sh
mpiexec -n 2 python tools/test_eh_deep_recovery_mpi.py BUILD_JSON NEW_OUTPUT --baseline --resume-midpoint --tensor
```
This reuses the common native recovery suite; it needs a source checkout and allocated resources.


## LD-only reference-frequency subpixel smoothing

Enable interface coefficient preparation explicitly:

```python
from salmonpy.subpixel import SubpixelSmoothing
from dataclasses import replace
smoothed = replace(calculation, smoothing=SubpixelSmoothing())
# PlaneWave: reference angular frequency comes from source.pulse.omega.
# Impulse has no carrier: specify SubpixelSmoothing(reference_omega=...).
```

This option is for native E/H FDTD (1D-like, 2D-like and 3D geometries) with
explicit LD materials. It never smooths SBE/TDDFT or user A-to-j constitutive
states. Interfaces must touch an LD phase; a constant dielectric may be its
other side. Without LD, or with unsupported magnetic/conductive background,
warn and retain the ordinary coefficients. Native response-point neighborhoods
are excluded from interface preparation. This does not smooth atomistic density
or change electronic current, its normalization or the Maxwell coupling protocol.

A dedicated Fortran module computes local planar fractions/normals for boxes
and ellipsoids and prepares passive tensor coefficients at the input carrier.
It preserves the original pole frequencies/damping, transforming each real
background and pole residue by frozen complex interface-field correctors.
The complex planar effective tensor matches at reference_omega, including normal
harmonic and tangential arithmetic response. No complex constant epsilon or
negative pole strength is introduced. No curl, ADE recurrence or E/J/P/Q/H update
order is replaced. The existing tensor solver handles the spatial coefficients.

This is a frequency-local dispersive approximation, NOT a claim of Meep's exact
implementation or second-order accuracy across an entire spectrum. Inspect
broadband response and grid/time convergence. Corners, distinct multiple
interfaces and singular reference mixtures retain local staircase coefficients
with warnings. LD interfaces with other LD materials are allowed independently
of their IDs. Source geometry overlap still follows last-shape precedence.
Vacuum PML/PEC treatment is retained; keep material away from PML.

The profile is salmonpy-subpixel.dat (version2), staged automatically by run()
and NativeEHField, bound to runtime/model/checkpoint identity. Export input_assets()
with input_text() when preparing a batch job. A matching audited native runtime
is required. Missing/corrupt profiles are errors, not silent fallback. Summary
salmonpy-subpixel-summary.json reports reference, global oscillator slots,
excluded response points, mixed/fallback sample counts, matching error and a
conservative background time-step bound. Exceeding this bound warns and continues;
finite evolution still requires its own checks. If a summary cannot be written,
warn and keep calculating.

Allocation scales with all active material pole slots times local grid volume,
including cell-owned tensors; poles are not currently merged by equal rates.
Inspect the native allocation estimate before using many distinct LD models.
Direct checkpoints include these coefficients and reject incompatible models;
MPI field repartitioning is not added by this option.

```sh
# Preview only: no calculation or output directory is started
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --subpixel --angle-degrees 30
# Explicit calculation inside allocated resources
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --subpixel --angle-degrees 30 --run
```

The sample's zero-strength LD background is intentional: it tests interface
geometry and polarization independently of fitted-material accuracy. Use an
actual dispersive LD model for a dispersive experiment. Keep smoothing=None
(default) for the unchanged staircase path.
