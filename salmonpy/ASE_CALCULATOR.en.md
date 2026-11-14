# ASE Calculator: energy, native forces and warm start

`salmonpy.ase_calculator.SalmonCalculator` is an optional independent module.
It supports energy (eV) and forces (eV/Å) for neutral, periodic, unpolarized
fixed-ion SCF: PZ, Libxc-PZ and Libxc-PBE, scalar KY or NC-UPF potentials.
Libxc requires a matching audited runtime. Check the potential's element/XC
provenance; no functional fallback is performed. TB-mBJ forces, stress,
free_energy, charged cells and partial periodic boundaries are unsupported.

## NLCC UPF forces

**Existing or unvalidated runtimes cannot provide Calculator NLCC forces.**
An ordinary `SalmonExecutable(path, sha256)` does not enable them. Use matching
PZ + KY potentials when optimizing without a validated NLCC runtime, and check
grid, k-point and SCF convergence.

Load a runtime with completed implementation validation:

```python
from salmonpy.runtime import load_runtime, load_validated_build

engine = load_runtime('my-runtime', required=['nlcc-forces'])
# Or use your trusted audited build and its completed validation record:
# engine = load_validated_build('work/build.json', 'checks/validation.json')
print(engine.nlcc_force_xc)  # ('pz', 'libxc_pz', 'libxc_pbe')
```

Pass this engine to Calculator to obtain NLCC forces in the validated scope:
CPU, unpolarized scalar NC-UPF, PZ/Libxc-PZ/PBE. The record binds executable SHA256,
source patch, XC/core finite-variation tests and a fixed Si8 application regression.
Bundling at a different location reruns the application check. Generate records
with trusted explicit build/validation commands; a feature label alone is insufficient.

**Implementation regression does not guarantee material force convergence.**
Check radial interpolation/cutoff boundaries, real-space grid, k points and SCF
convergence for your system. The compatibility patch corrects the old Libxc-NLCC
energy weighting/GGA gradients and adds the explicit core force. Existing binaries
are not updated automatically.

Without matching validation, a truthy/unrecognized `core_correction` flag or any
`PP_NLCC` element disables forces before SCF with `PropertyNotImplementedError`.
Energy remains available without caching forces. Availability alone does not
certify energy accuracy. Lower-level raw `GroundStateResult` force access also
does not certify its runtime.

```python
from salmonpy.ase_calculator import SalmonCalculator
from salmonpy.dft import SCFSettings

atoms.calc = SalmonCalculator(
    executable=engine, directory='new-calculations',
    pseudopotentials=[pseudo], states=32, grid=(16, 16, 16), kmesh=(2, 2, 2),
    xc='pz', scf=SCFSettings(max_iterations=500, tolerance=1e-11),
    threads=1, warm_start=True,
)
energy_ev = atoms.get_potential_energy()
forces_ev_per_angstrom = atoms.get_forces()
```

Prepare `atoms`, a matching PZ KY `pseudo`, and the audited `engine` first.
These grid/band/k-point values are small Si connection-test values, not general
accuracy defaults. `electrons` may be omitted: the neutral pseudopotential valence
sum is inferred and explicit counts must agree. States/grid/kmesh remain required.
See [ASE setup and optical response](ASE.en.md).

## Wavefunction reuse

The default `warm_start=False` preserves fresh-SCF behavior. With True, a new
geometry reuses wavefunctions only from the last successful compatible result.
Cell, grid, k points, ordered species, potential bytes/local channels, electron
count, orbital count and XC must match. The executable is fixed on the Calculator.
Incompatible model conditions cause a fresh SCF; corrupted saved files cause an
error rather than silently falling back. Failure clears cached results and the
previous initial-guess candidate. The old source directory is never modified.

Each calculation owns a new output directory. Wavefunctions and occupations are
copied into its initial_guess directory, density is recalculated for the new
geometry, and SCF iteration count is reset. Old ionic coordinates and mixing
history are not reused. Source identity is recorded in ground-state.json.
The ASE-independent equivalent is
`model.run(engine, NEW_DIRECTORY, initial_guess=previous_ground_state)`;
this explicit path rejects incompatible guesses before launching native code.

Warm start is not guaranteed to reduce iterations or reach the same local minimum.
Compare fresh and warm energies/forces at the same SCF tolerance for your system.
Check grid convergence and rigid translations: real-space egg-box effects mean
coarse-grid finite-difference energy derivatives are not a force reference.

## Scope and lifetime

Repeated queries for an unchanged structure use ASE's cache. Changed structures
run SCF again. Settings are immutable; create a new Calculator to change them.
`last_ground_state` is the last successful saved GS, and is cleared on failure.
Use a separate Calculator per concurrent calculation. Directory retention means
long optimizations can use substantial disk space.

Fixed-cell ASE optimization/Born–Oppenheimer MD is distinct from coupled
Ehrenfest electron/ion propagation. Isolated systems, non-collinear spin and
Ehrenfest remain first-release targets requiring native interface work.
Crystalline TDDFT uses **transverse** (`trans_longi='tr'`) as its standard;
Longitudinal is an explicitly chosen field component, not a default.
