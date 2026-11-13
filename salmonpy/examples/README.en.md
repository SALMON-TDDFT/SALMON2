# Python samples

[日本語](README.md) · Topic: `examples` · Paired revision: `2026-11-13.1`

Use the same Python environment in which the wheel was installed. All samples
below ship in the wheel. The code is shared by both language editions.

| Module under `salmonpy.examples` | Purpose / additional input |
|---|---|
| `ld_response` | LD time response and permittivity; bundled native kernel |
| `pulse_spectrum` | Pulses and Fourier spectrum; no solver |
| `geometry_preview` | Material IDs, middle slice, PML contact; new NPZ path |
| `sbe_response` | Synthetic two-band A-to-J and independent states; new directory |
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
python -m salmonpy.examples.outputs_roundtrip new-ld-output
```

These five examples need no external input dataset. Synthetic SBE data illustrates
the API, not an actual material. Core units are atomic units; A means A/c and
crystalline electronic charge current is −J. LD is E-driven.

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
checkout. BLAS threads also count against your CPU budget. Python solver launches
use one MPI rank. Sample grid/time values are demonstrations, not convergence or
memory guarantees. TDDFT's short pulse is a connection smoke test. The Calculator
sample does not automatically start an optimization. Adapt plotting assumptions
for datasets other than the specified cubic Si example.

Published SALMON inputs can be found at
[SALMON-inputs](https://github.com/SALMON-TDDFT/SALMON-inputs).
Check the exact input revision, solver version, required files and wrapper scope
before adapting one; these are not certified SALMONpy fixtures.

Public SALMON VR export: [salmon_vr_export.py](salmon_vr_export.py), see [guide](../SALMON_VR.en.md).
