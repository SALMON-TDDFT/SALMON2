# Python samples

[日本語](README.md) · Topic: `examples` · Paired revision: `2026-11-13.5`

Use the same Python environment in which the wheel was installed. All samples
below ship in the wheel. The code is shared by both language editions.

LD and built-in SBE use the native library bundled in the wheel. Pulse creation and analysis need no solver.
For SALMON input exported by `FieldSamples.write_salmon`, use `ae_shape1='input'` and `unit_system='au'`.

| Module under `salmonpy.examples` | Purpose / additional input |
|---|---|
| `ld_response` | LD time response and permittivity; bundled native kernel |
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

These five examples need no external input dataset. Synthetic SBE data illustrates
the API, not an actual material. Core units are atomic units; A means A/c and
crystalline outputs are the electric (charge) current density j (SALMON matter current = −j). LD is E-driven.

```sh
python -m salmonpy.examples.dft_si EXECUTABLE SHA256 SI_RPS NEW_GS --threads 1
python -m salmonpy.examples.tddft_response EXECUTABLE SHA256 NEW_GS NEW_RT --threads 1
python -m salmonpy.examples.bands_from_gs EXECUTABLE SHA256 NEW_GS NEW_BANDS --threads 1
python -m salmonpy.examples.ase_bands EXECUTABLE SHA256 NEW_GS NEW_ASE_BANDS --threads 1
python -m salmonpy.examples.ase_energy_forces EXECUTABLE SHA256 SI_RPS NEW_CALCULATIONS --threads 1
python -m salmonpy.examples.fdtd_dielectric EXECUTABLE SHA256 NEW_FDTD
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

## Couple native E/H FDTD and internally parallel TDDFT

[mpi_coupled_tddft.py](mpi_coupled_tddft.py) uses two field ranks plus two two-rank
material groups. Follow the [host-library setup and time/state contract](../COUPLED.en.md).

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 6 python -m salmonpy.examples.mpi_coupled_tddft BUILD_JSON GS_BUNDLE RT_INPUT NEW_OUTPUT
```

`response.json` stores two-point A/E/J and both clocks; the field directory contains
the manifest and observations. One native lifetime per process; another trial needs a fresh MPI job.

For practical TDDFT, prefer 1D propagation/layered structures; recommend 3D electronic feedback with SBE or a user-written current solver. This 40³/two-point example is a connection regression, not a 1D reduction or recommended production 3D TDDFT setup.
