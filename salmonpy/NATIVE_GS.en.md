# Retained standard ground states

[日本語](NATIVE_GS.md) · Topic: `retained_native_gs`

salmonpy.native_gs keeps several standard DFT ground states in one process.
Each owns its input, orbitals, density, occupations, eigenvalues, forces, XC and
MPI children. Creating a different cell/grid/k/electron/atom/XC calculation leaves
previous owners readable. Ordinary DFT.run retains its separate-process behavior.

```python
from types import SimpleNamespace
from mpi4py import MPI
from salmonpy.dft import SCFSettings
from salmonpy.native_gs import NativeGSBackend

# dft_a/dft_b are immutable DFT objects; use a work root shared by subgroup ranks.
context = SimpleNamespace(solver_comm=MPI.COMM_WORLD, closed=False, _dependents=0)
backend = NativeGSBackend(build_record, work_root,
                          decomposition=(1, MPI.COMM_WORLD.size, 1, 1, 1))
with backend.open(context) as session:
    a = session.run(dft_a, scf=SCFSettings(), material_id="film-a")
    b = session.run(dft_b, scf=SCFSettings(), material_id="film-b")
    energy_a = a.energy_hartree
    density_a = a.array("density")
    eigen_b = b.array("eigenvalues")
    a.close()
    energy_b = b.energy_hartree
# All states close before native detach; MPI_Finalize remains the host's job.
```

Decomposition is (k, orbital, x, y, z); its product equals subgroup size. A host
may pass its existing split communicator as context.solver_comm. All ranks call
matching operations, inputs, material IDs and handles in order on the main thread.
Lightweight LD/SBE components remain individually usable.

session.run computes SCF. state.info returns total Hartree energy, residual,
tolerance, iteration count, the native recorded convergence flag, atom and band
counts. Unconverged SCF warns and retains the last iterate. material_id and the
immutable calculation remain attached to the state.

state.array(quantity) returns NativeGSArray(values, lower_indices, quantity, unit).

| quantity | Meaning | Unit |
|---|---|---|
| density | Total valence-electron density | bohr^-3 |
| occupations | Occupations | dimensionless |
| eigenvalues | Single-particle energies | Hartree |
| forces | Atomic forces | Hartree/Bohr |
| orbitals | Real/complex orbitals returned as complex arrays | native orbital amplitude |

values is a read-only copy, independent of native storage and valid after close.
Axes and lower_indices preserve native Fortran layout. Density/orbitals are local
to each rank; orbitals may include halos. Orbital axes are
(x,y,z,spin,orbital,k,image); eigen/occupation axes are (orbital,k,spin), and forces
are (component,atom). The latter three arrays are natively replicated. There is
no automatic global concatenation; do not double-count replicas.

Retained owners and requested copies consume memory. Read large orbital arrays
only when needed. session.resources reports states, live child communicators and
management slot high-water. Closing one owner preserves siblings; closing all
returns live children to zero. Retained management slots are distinct from leaks.

Python scope matches existing DFT settings: periodic full-k unpolarized fixed-ion
PZ/Libxc-PZ/PBE. Fast-PR HSE, SOC/+U/DC GS, direct in-memory GS→TDDFT handoff,
checkpoint/restart, concurrent threads and simultaneous GS/TDDFT are separate.
One native electronic attachment lifetime per process. Info, arrays and close are
also collective. Close every owner before host MPI finalization.

Use an explicitly built/audited standard runtime with sp_gs_protocol()==1. Older
prebuilt runtimes may lack it. Imports never acquire or compile runtimes.
tools/build_embedded_tddft.py builds from audited sources/application objects.
Darwin linking keeps NumPy BLAS symbols separate while sharing MPI COMMON sentinels;
Fortran host probes use mpi_link_command for the same convention. Linux linking
stays unchanged. Current measurements are Apple Silicon/GNU/OpenMPI; other
platforms require their own verification.

Work directory/input/pseudopotential visibility and content are checked on every
rank before native calls. General Fortran STOP/allocation/MPI/process failure
recovery is not promised. Numerical update order and warning-only SCF policy stay.
The example `mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD NEW_OUTPUT`
uses coarse grids/four iterations for ownership verification, not converged
material properties. It writes normal native diagnostics/progress, but no automatic
GS restart bundle. See [the native contract](native_adapters/standard_gs.h).
