# ASE structures, bands and optical response

The optional `salmonpy.ase` module copies ASE structures and paths into the
ASE-independent DFT/band/TDDFT interfaces. Importing it does not import ASE;
calling ASE operations requires the optional dependencies:

```sh
python -m pip install "ase>=3.26,<4" "f90nml>=1.4,<2"
```

Install a compatible SALMONpy wheel first. Native DFT/TDDFT/bands also need an
explicitly installed, audited external runtime (see [runtime setup](RUNTIME.en.md)).
Nothing is compiled or downloaded on import. The bridge does not change native
solver order. Ordinary LD/SBE/TDDFT use does not require ASE.

## Neutral structures and electron counts

`dft_from_atoms(atoms, pseudopotentials=..., states=..., grid=..., kmesh=...)`
returns an independent DFT model without running a calculation. Cell rows and
positions are converted from Å using SALMON's 0.52917721067 Å/bohr.
It preserves the cell and ordered atoms; no primitive-cell conversion is implicit.

`electrons=None` uses the sum of pseudopotential valence charges over atoms.
KY uses the fourth field in its first header record; NC-UPF uses
`PP_HEADER/z_valence`. `Pseudopotential.valence_electrons` exposes this value.
An explicit count must match the neutral sum. Invalid/noninteger charges are
rejected, never rounded or inferred from atomic numbers. The supplied atomic
number and potential/XC identity still need verification; KY does not identify
its element in this charge header. Odd-electron closed shells are unsupported.

`states`, `grid` and `kmesh` remain required. No automatic setting can guarantee
convergence. Choose initial values against the purpose and memory/CPU budget,
then check the target observable against grid, k-point, band and time-step changes.
The ASE Calculator also accepts omitted electrons, resolving/checking them for
each structure before native execution. It supports PZ/Libxc-PZ/Libxc-PBE with scalar KY/NC-UPF potentials, energy and forces only.
Libxc needs a matching runtime. Validate potential/XC identity. Stress is not supported;
each geometry retains a separate GS directory. See [Calculator and warm start](ASE_CALCULATOR.en.md).

## Atoms to A-driven current

```python
from ase.build import bulk
import numpy as np
from salmonpy.ase import BOHR_ANGSTROM, tddft_from_atoms
from salmonpy.dft import Pseudopotential, SCFSettings
from salmonpy.runtime import load_runtime
from salmonpy.fields import AcosPulse
from salmonpy.crystal import respond

engine = load_runtime('installed-runtime', required=('dft', 'tddft', 'pz', 'single-rank'))
atoms = bulk('Si', 'diamond', a=10.26*BOHR_ANGSTROM, cubic=True)
pseudo = Pseudopotential.from_file('Si_rps.dat', atomic_number=14, local_channel=2)
drive = AcosPulse(amplitude=.001, omega=.2, duration=.4,
                 polarization=(0, 0, 1)).sample(np.arange(21)*.02)
with tddft_from_atoms(atoms, 'new-gs', executable=engine,
                      pseudopotentials=[pseudo], states=32,
                      grid=(16, 16, 16), kmesh=(2, 2, 2),
                      scf=SCFSettings(max_iterations=500, tolerance=1e-10),
                      max_steps=20, threads=1) as initial:
    response = respond(initial, drive, output_directory='new-current')
    with initial.new_state(dt=.02) as state:
        j = state.step([0., 0., 1e-4])
```

This is a short connection example, not converged optical physics. For verified
starter inputs/potential/reference use the [Si tutorial](SI_INTRO.en.md).
The complete packaged example is:

```sh
python -m salmonpy.examples.ase_tddft RUNTIME_DIRECTORY NEW_OUTPUT_DIRECTORY --verify-independence
```

It runs 8-atom Si, 16³ grid, 2³ k points and 20 steps with one thread, writes
`gs/`, `rt/`, `current.npz` (time, A, current) and `report.json`, and compares against
the starter reference with tolerances. The optional independence check runs an
opposite waveform then replays the original from the same initial state.
The example does not measure or enforce peak memory; resource status is separate
from calculation/validation status.

`tddft_from_atoms` returns the existing context-managed `TDDFTInitialState`.
Each `new_state` or `respond` starts an independent history. Closing the initial
snapshot leaves the saved GS and logs intact; reuse via `GroundStateResult.open`.
A is A/c in atomic units. Returned current is **charge-current density j**,
minus SALMON's electronic matter current J; do not negate again.

Scope: neutral, periodic, unpolarized, fixed ions, one MPI rank;
PZ/Libxc-PZ/Libxc-PBE with a matching runtime and pseudopotential. This entry does
not add TB-mBJ GS creation, isolated systems, non-collinear spin or ionic motion.

## Band and geometry bridge

`atoms_from_initial(initial)` reads the original GS geometry (inline or relative
coordinate file; Cartesian or fractional; au or A_eV_fs) into new Atoms.
`run_bands(initial, atoms, path, directory, reference_ev=..., **run_options)`
validates structure and passes ASE BandPath points unchanged to the native band
runner. Ordered species must match, positions are equivalent modulo lattice
translations, and cell/position tolerance is 1e-8 Å. No provenance of binary GS
files is inferred beyond the supplied input and manifest.

The returned `ASEBandRun` retains the native result and ASE `band_structure`.
It saves `ase-bands.json` and `ase-input.json`; `.plot(...)` delegates to ASE.
For conversion only, use `explicit_path` or
`to_ase_band_structure(bands, path, reference_ev=...)`. The energy reference is
explicit, never an inferred Fermi level. Arrays change from (k, band, spin) Hartree
to (spin, k, band) eV using 27.211386245988 eV/Hartree. Actual Cartesian k points
must match the supplied path/cell. Headless plotting needs an Agg backend.

The standard crystalline TDDFT configuration is **transverse** (`trans_longi='tr'`).
Longitudinal is an explicit alternative field component, not the default.
Isolated systems, non-collinear spin and Ehrenfest remain first-release targets;
they are not yet implemented by this ASE entry.
