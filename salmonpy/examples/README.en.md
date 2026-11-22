# Python samples

[日本語](README.md)

[Si starter: wheel to current](../SI_INTRO.en.md) includes a fixed pseudopotential, input and reference, and explicitly selects a prebuilt runtime.

Use the same Python environment in which the wheel was installed. All samples
below ship in the wheel. The code is shared by both language editions.

LD and built-in SBE use the Fortran/C library bundled in the wheel. Pulse creation and analysis need no solver.
For SALMON input exported by `FieldSamples.write_salmon`, use `ae_shape1='input'` and `unit_system='au'`.

## Start with these three

1. **LD response**: `python -m salmonpy.examples.ld_response`. Try material response using only the wheel.
2. **SBE response**: `python -m salmonpy.examples.sbe_response new-sbe-output`. Learn A-to-j and independent states using a synthetic model, not real-material optical data.
3. **Si current**: follow the [Si starter](../SI_INTRO.en.md) with a matching prebuilt runtime. This short wheel-to-current check does not certify physical convergence.

These three are enough for a first visit. The table below indexes further tasks; MPI coupling examples are advanced. Use new output locations outside the source checkout.

| Module under `salmonpy.examples` | Purpose / additional input |
|---|---|
| `custom_checkpoint` | User-calculation implementation capability/identity/compatibility checks and independent continuation; new output directory, no external solver |
| `ld_response` | LD time response and permittivity; bundled compiled calculation routine |
| `pulse_spectrum` | Pulses and Fourier spectrum; no solver |
| `geometry_preview` | Material IDs, middle slice, PML contact; new NPZ path |
| `sbe_response` | Synthetic two-band A-to-J and independent states; new directory |
| `field_coupling` | Field components (tr / lo / 1D propagation) coupled to synthetic SBE; new directory |
| `sbe_from_salmon` | Real-material loader; `_k/_eigen/_tm.data` and column lattice |
| `outputs_roundtrip` | Save/read selected LD arrays with units; new directory |
| `fdtd_dielectric` | Dielectric propagation; audited external SALMON |
| `dft_si` | Eight-atom cubic Si GS and forces; SALMON and four-valence-electron PZ/KY Si potential |
| `tddft_response` | Two sequential independent responses from one saved GS |
| `bands_from_gs` | Cubic Si GS to bands; SALMON and Matplotlib |
| `plot_bands` | Existing 17-point/24-band/32-electron cubic Si output; Matplotlib |
| `ase_bands` | Cubic Si GS and ASE path; ASE, f90nml, Matplotlib |
| `ase_energy_forces` | PZ/KY energy and forces; ASE and full SALMON |

```sh
python -m salmonpy.examples.ld_response
python -m salmonpy.examples.pulse_spectrum
python -m salmonpy.examples.geometry_preview new-grid.npz
python -m salmonpy.examples.sbe_response new-sbe-output
python -m salmonpy.examples.field_coupling new-field-coupling-output
python -m salmonpy.examples.outputs_roundtrip new-ld-output
```

These six examples need no external input dataset. Synthetic SBE data illustrates
the API, not an actual material. Core units are atomic units; A means A/c and
crystalline outputs are the electric (charge) current density j (SALMON matter current = −j). LD is E-driven.

```sh
python -m salmonpy.examples.dft_si EXECUTABLE SHA256 SI_RPS NEW_GS --threads 1
python -m salmonpy.examples.tddft_response EXECUTABLE SHA256 NEW_GS NEW_RT --threads 1
python -m salmonpy.examples.bands_from_gs EXECUTABLE SHA256 NEW_GS NEW_BANDS --threads 1
python -m salmonpy.examples.ase_bands EXECUTABLE SHA256 NEW_GS NEW_ASE_BANDS --threads 1
python -m salmonpy.examples.ase_energy_forces EXECUTABLE SHA256 SI_RPS NEW_CALCULATIONS --threads 1
python -m salmonpy.examples.fdtd_dielectric RUNTIME NEW_FDTD
python -m salmonpy.examples.sbe_from_salmon BUNDLE SYSNAME NEW_SBE --bands 32 --kpoints 8 --electrons 32 --lattice lattice-columns-bohr.txt
```

Replace every uppercase placeholder. Paths containing spaces need quotes. SBE
counts and lattice must match the actual dataset. Argument-based examples support
`--help` (some descriptions remain Japanese). No sample downloads potentials.

Use a compatible audited full SALMON, not an arbitrary installed binary. See
[setup and scope](../README.en.md). Choose a new output location outside the code
checkout. BLAS threads also count against your CPU budget. Python DFT/bands and process-isolated `salmonpy.tddft` solver launches
use one MPI rank. Internally parallel TDDFT uses the separate `mpi_coupled_tddft` entry.
The MPI FDTD example explicitly launches its spatial ranks. Sample grid/time values are demonstrations, not convergence or
memory guarantees. TDDFT's short pulse is a connection smoke test. The Calculator
sample does not automatically start an optimization. Adapt plotting assumptions
for datasets other than the specified cubic Si example.

Published SALMON inputs can be found at
[SALMON-inputs](https://github.com/SALMON-TDDFT/SALMON-inputs).
Check the exact input revision, solver version, required files and wrapper scope
before adapting one; these are not certified SALMONpy fixtures.

Public SALMON VR export: [salmon_vr_export.py](salmon_vr_export.py), see [guide](../SALMON_VR.en.md).

## Multiple materials and MPI

Replace uppercase arguments with local paths. `RUNTIME` is a managed **serial** runtime;
`MPI_BUILD_JSON` is the completed MPI build's `build.json`. They are not interchangeable.
`GS_BUNDLE` is a self-contained directory with RT input, pseudopotentials and restart data;
`RT_INPUT` is its input filename. Match the GS functional, grid and other settings.

| Sample | Allocation | Output |
|---|---|---|
| [fdtd_layers.py](fdtd_layers.py) | Serial FDTD, two synthetic LD layers, one thread | Observation logs, input, run manifest |
| [fdtd_mpi_layers.py](fdtd_mpi_layers.py) | Four field ranks, one thread each; three LD materials, 768 matter cells | Observation logs, input, manifest with MPI decomposition |
| [mpi_response.py](mpi_response.py) | One field rank + two two-rank material groups = five ranks; one material/point per group | JSON A/c and charge-current density for point IDs 10/20 |
| [mpi_tddft.py](mpi_tddft.py) | One field rank + two one-rank material groups = three ranks; two serial TDDFT workers | JSON A/c and charge-current density for point IDs 1/2 |

```sh
python -m salmonpy.examples.fdtd_layers RUNTIME NEW_LD_DIRECTORY
python -m salmonpy.examples.fdtd_mpi_layers MPI_BUILD_JSON NEW_MPI_DIRECTORY --grid 2 2 1 --threads 1
mpiexec -n 5 python -m salmonpy.examples.mpi_response NEW_CUSTOM_JSON
mpiexec -n 3 python -m salmonpy.examples.mpi_tddft RUNTIME GS_BUNDLE RT_INPUT NEW_TDDFT_JSON
```

Start FDTD samples from **fresh serial Python**. The response samples start Python under MPI.
Replace `MyMaterial.prepare` with your loader and `MyResponse.create_state` with your response
implementation; group collectives use the supplied `solver_comm`. Both TDDFT points start
from the same GS but keep independent excitation histories. Memory limits are not checked
automatically; budget for all states, worker processes and BLAS threads.

Inspect `salmonpy-run.json` for FDTD completion, rank count, decomposition and executable.
Response JSON `history` contains atomic-unit `time`, `A_over_c`, `J_charge`; rows follow
`point_ids`. It is postprocessing data, not a wavefunction/full-state restart file.

LD samples use finite synthetic structures; the custom response uses a manufactured memory
law. The TDDFT sample is a short prescribed-A connection test. None is a converged multilayer
optics calculation. The response samples do not propagate FDTD or provide inner-MPI TDDFT.
See [material interfaces](../FDTD_MATERIALS.en.md) and [distributed scope](../DISTRIBUTED.en.md).

## Couple FDTD called directly from Python and internally parallel TDDFT

[mpi_coupled_tddft.py](mpi_coupled_tddft.py) uses two field ranks plus two two-rank
material groups. Follow the [calling program-library setup and time/state contract](../COUPLED.en.md).

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 6 python -m salmonpy.examples.mpi_coupled_tddft BUILD_JSON GS_BUNDLE RT_INPUT NEW_OUTPUT
```

`response.json` stores two-point A/E/J and both clocks; the field directory contains
the manifest and observations. One solver initialization per process; another trial needs a fresh MPI job.

For practical TDDFT, prefer 1D propagation/layered structures; recommend 3D electronic feedback with SBE or a user-written current solver. This 40³/two-point example is a connection regression, not a 1D reduction or recommended production 3D TDDFT setup.

## ASE to Si current

`python -m salmonpy.examples.ase_tddft RUNTIME_DIRECTORY NEW_OUTPUT_DIRECTORY --verify-independence`

Requires ASE/f90nml and an audited DFT/TDDFT/PZ runtime. See the [ASE guide](../ASE.en.md). Infers neutral electron count and runs a short Si connection test with one thread, not converged optical physics.

## Snapshot replay

```sh
python -m salmonpy.examples.snapshot_restore
```

Checks LD/SBE and all four Python fields; installed calculation routines, no full SALMON required. See [snapshots](../SNAPSHOTS.en.md).

Custom A-driven models can use one serial/MPI contract and an explicit conformance checker: [authoring guide](../MATERIALS.en.md), `python -m salmonpy.examples.custom_material`.

Validate material IDs and custom settings with [material_config.py](material_config.py). No extra data; configuration only, no simulation.

See [CONFIG.en.md](../CONFIG.en.md) / [CONFIG.md](../CONFIG.md).

Stage local input files with [prepare_inputs.py](prepare_inputs.py); see the [guide](../PREPARATION.en.md). Supply a synthetic coefficient JSON and a new output directory.

Prepare a structure/grid input bundle with [prepare_structure.py](prepare_structure.py); see [guide](../STRUCTURE.en.md). Requires ASE.

Run two LD materials with explicit IDs: [configured_ld.py](configured_ld.py). Requires a new output directory and installed Fortran/C library.

Run configured Si DFT→TDDFT: [configured_si.py](configured_si.py), [guide](../DFT_CONFIG.en.md). Requires an existing runtime and ASE/f90nml.

## Additional connection and input samples

| Sample | Purpose | Requirements / arguments |
|---|---|---|
| [field_inputs.py](field_inputs.py) | A/E/H input contract and custom response | LD calculation routine |
| [point_ld.py](point_ld.py) | LD with material and point identities | LD calculation routine; new output |
| [point_tddft_restart.py](point_tddft_restart.py) | Per-point TDDFT restart | Catalog, runtime, waveform; new output |
| [catalog_tddft.py](catalog_tddft.py) | TDDFT from a catalog | Catalog, runtime, waveform; see --help |
| [mpi_configured_native.py](mpi_configured_native.py) | Configured LD/TDDFT MPI response | 3 ranks; catalog, waveform; TDDFT runtime |
| [mpi_recorded_material.py](mpi_recorded_material.py) | Recorded custom/SBE serial–MPI comparison | 3 ranks; optional SBE catalog |
| [mpi_point_identity.py](mpi_point_identity.py) | Point and material identities over MPI | 3 ranks; new output |
| [mpi_multiscale_tddft.py](mpi_multiscale_tddft.py) | 1D field coupled to independent TDDFT groups | 3 ranks with --material-ranks 1; solver build and GS |
| [isolated_response.py](isolated_response.py) | Isolated E-driven dipole response | Validated E-worker, validation record, GS bundle, E NPZ |
| [dipole_radiation.py](dipole_radiation.py) | Dipole-to-current mapping for FDTD | No data, or --archive and --cell-volume; no FDTD propagation |

MPI samples also display `--help` from ordinary Python. Computation requires mpi4py and compatible MPI. Check rank options before launching; some defaults exceed a four-core budget.

[mpi_field_inputs.py](mpi_field_inputs.py): three ranks with distinct A/E/H contracts and synthetic currents; verifies rejection before updates.

[field_sampling.py](field_sampling.py): explicit spatial weights and staggered-time interpolation of synthetic fields; no solver/data required.

[mpi_packet_coupling.py](mpi_packet_coupling.py): synthetic packet feedback with one-dimensional fields or FDTD called from Python; three ranks (solver also supports four).

- `molecular_field_selection.py`: explicit own-source vacuum field subtraction on synthetic samples; no solver required.

Larger parallel checks: `mpi_mixed_native` supports `--field-shape`. See [COUPLED.en.md](../COUPLED.en.md) for5–15-rank layouts and measurement scope. The4GiB/4-core beginner configuration remains available.

Offline convergence planning: `python -m salmonpy.examples.convergence_plan NEW_JSON`; see CONVERGENCE.en.md. Explicit illustrative PP radius only; no material computation.

## Start FDTD with these three examples

1. `fdtd_dielectric RUNTIME NEW_OUTPUT`: open a verified bundle through `load_runtime` and check dielectric propagation.
2. `fdtd_layers RUNTIME NEW_OUTPUT`: try two LD layers and material-ID-independent interfaces.
3. `fdtd_mpi_layers MPI_BUILD_JSON NEW_OUTPUT`: use an audited MPI-capable build record for spatial parallelism.

`RUNTIME` is a directory containing `runtime.json`, not a bare executable. Existing applications using `load_runtime` keep that entry unchanged. The MPI example is advanced and requires its separate capability-bearing build record. Start from fresh serial Python and use new output directories.

`fdtd_periodic_slab RUNTIME NEW_OUTPUT [--ld]` demonstrates a grid-aligned slab with transverse periodic boundaries. See [FDTD materials](../FDTD_MATERIALS.en.md) for the C-reference Fresnel test.

### Serial TB-mBJ ground states

Use `DFT(..., xc="tbmbj", cval=None)` or `dft_from_atoms(..., xc="tbmbj", cval=1.2)`. Material JSON recipes accept xc and optional cval. None selects the automatic coefficient; a finite value greater than 1e-10 selects a fixed coefficient. The coefficient is retained in GS records, warm-start compatibility checks and TDDFT inputs. TB-mBJ is a model potential: total energies and forces are unvalidated diagnostics, with a warning rather than a run prohibition. ASE Calculator optimization and solver parallel SCF support are unchanged. Run `python -m salmonpy.examples.dft_tbmbj --help` for the serial example. Validation here covers synthetic solver XC evaluation and interface contracts, not real-material accuracy.

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

Constrained LD fit: `python -m salmonpy.examples.ld_fit --poles 1 2 3`. Compare residuals, acceptance, timestep advice and provenance for explicitly selected pole counts on synthetic data.

Offline measured optical data: `python -m salmonpy.examples.optical_materials`. Lists pinned CC0 measurements and loads stored coefficients without downloading/refitting.

`python -m salmonpy.examples.ld_response --anisotropic`: synthetic reciprocal tensor material response, not solver anisotropic FDTD propagation.

`fdtd_waveplate`: quarter/half-wave synthetic FDTD; default preview only, --run starts computation. See FDTD_MATERIALS.en.md.

For optional LD-only interface smoothing, add `--subpixel` to `fdtd_waveplate`. It matches the carrier locally; it does not transform SBE/TDDFT. Default is off. Preview first and explicitly pass `--run` inside allocated resources.

`optical_postprocess` compares duration/windows and weighted k contributions on stored arrays. Preview is default; `--run` performs postprocessing only, never launches a solver. See OPTICAL_TRAJECTORIES.en.md.

`mpi_runtime_receipt` measures explicit within-job validation I/O. Preview is default; even --run never launches SALMON. See RUNTIME.en.md.

Mixed MPI sessions route A/field packets only to owning points; J return order is unchanged. See communication_report for the ownership layout.

With a matching spatial-current solver runtime, mpi_mixed_native also routes J by field ownership. Only owned coupling-J rows are meaningful off the field leader; use leader output for full-point analysis. See COUPLED.en.md.

Preparation only: `python -m salmonpy.examples.prepare_fdtd NEW_DIRECTORY --mpi-grid 2 1 1` writes FDTD inputs and stops, without solver execution or submission. See PREPARATION.en.md.

Electronic input preparation only: `python -m salmonpy.examples.prepare_electrons dft Si_rps.dat NEW_DIRECTORY` or `tddft rt.inp NEW_DIRECTORY --asset pseudo1_rps.dat Si_rps.dat`. No solver, GS discovery or submission. See PREPARATION.en.md.

`fdtd_energy`: diagnose macro-point work, R/T and energy balance from stored collocated E/j and separated E/H ports; no solver launch.

`fdtd_periodic_slab --native-planes`: calculation object-local plane streams, automatic colocation and integrated flux; unmatched runtime warns without discarding the calculation.

The mixed solver sample optionally accepts `--planes-x FRONT BACK --plane-stride 16`: transverse-periodic planes captured only after successfully completed steps, saved once as plane-series. Disabled by default.

## Molecular feedback to FDTD

`molecular_feedback.py` is a complete MPI example reusing an existing GS. Required calculation programs and invocation are documented in [Molecular/FDTD coupling](../DIPOLE_RADIATION.en.md). It does not implicitly run SCF or a build.

`isolated_response.py`: `--stop-after 7 --direct-checkpoint NEW` saves the molecule at step7; a separate invocation with the same input can use `--resume-direct NEW` to evolve only remaining E samples. Both calls need the same capacity/runtime and a revision5 validation record.
