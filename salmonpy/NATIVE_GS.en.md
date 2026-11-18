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
PZ/Libxc-PZ/PBE. Fast-PR HSE, SOC/+U/DC GS,
full SCF checkpoint/restart and threaded GS/TDDFT are separate; sequential
GS/RT interleaving and portable array export are described below.
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

## Existing GS and in-memory TDDFT

`session.load(directory, input_file="rt.inp", material_id="silicon")` reads an
existing GS checkpoint using the supplied matching input, without SCF and without
writing a new GS restart bundle. Every rank must see unchanged source files.
Legacy binary integrity is checked; physical provenance is the caller's responsibility.
Imported `info.converged` is `None`, not a claim of converged SCF. Its residual and
tolerance fields are zero sentinels, not measured SCF convergence diagnostics.
Forces are not computed on import; requesting them raises an error.

```python
with backend.open(context) as session:
    gs = session.load(existing_bundle, material_id="silicon")
    model = gs.tddft_model(dt=0.003, max_steps=100)
    gs.close()
    first = model.new_state("first")
    second = model.new_state("second")
    j = first.step([0.0001, 0, 0])  # transverse, A/c -> charge current j
    first.checkpoint(checkpoint_path)  # each rank needs a distinct path
    # A fresh state's checkpoint(..., restore=True) restores without time replay.
    first.close()
    second.close()
    model.close()
```

The model owns a copy of initial orbitals, occupations and input values. Closing GS
preserves the model; new states are independent and start at the same initial GS.
A retained GS produced by `session.run` provides the same `tddft_model` method.
GS imports/solves and RT steps may be interleaved sequentially within one session.
Session close retires all RT models/states before releasing the native attachment.
Model creation currently inherits the GS subgroup and decomposition, supports the
standard periodic full-k unpolarized fixed-ion path, and rejects too many orbital
ranks for occupied RT bands. No automatic MPI repartitioning or thread concurrency.
Use a runtime exporting `sp_tddft_memory_gs_protocol()==1` and `sp_gs_import`.

```sh
mpiexec -n 2 python -m salmonpy.examples.memory_gs BUILD EXISTING_GS NEW_WORK
```

## Explicit portable export from a retained GS

For a regular-mesh state created by `session.run`, all subgroup ranks call:

```python
from salmonpy.gs_archive import GroundStateArchive

archive = gs.export_archive("new-portable-gs", orbital_block=16)
# Closing or subsequently changing native owners does not change the archive.
gs.close()
if archive is not None:  # subgroup root only; workers return None
    archive = GroundStateArchive.open("new-portable-gs")
    u = archive.wavefunctions(0)  # k=0, (band,z,y,x), read-only
```

This uses the existing GS/2 format and needs no new native ABI or compilation.
It preserves orbital/occupation/eigenvalue values, removes orbital halos using
physical density bounds, and reconstructs density in band order from occupations
and orbitals. The original mixed density is not exported. Material ID, GS identity,
runtime hash, recorded convergence, iteration count and pseudopotential hashes/local
channels are retained as provenance. Pseudopotential bytes are not embedded.
Unconverged SCF only warns and exports the last iterate. Existing destinations
are preserved; failed writes clean their newly claimed incomplete archive.

Current export scope is a DFT-owned, periodic unpolarized fixed-ion regular-k
state with a host MPI subgroup. Both Gamma and multiple k points are supported. Existing-file imported owners lack an immutable DFT model
and cannot use this entry yet; native GS-handle restoration remains future work.
These are export limits, not restrictions on normal GS/TDDFT calculations.
The first implementation is verified with synthetic real/complex-shaped native
arrays, asymmetric halos and noncubic grids; real retained-GS export and resulting
TDDFT parity have not yet been measured. No new SCF or real GS collection was
performed for this change.

This is an array interchange export, not a full SCF restart: mixing history,
forces, solver resources and original mixed density are not saved.
`GroundStateArchive.open` reads portable arrays; it does not restore a native GS
handle. Native import still uses the separately documented, explicitly built
paired converter and its runtime/potential restrictions in [GS_ARCHIVE.en.md](GS_ARCHIVE.en.md).
The optional sample `python -m salmonpy.examples.retained_gs BUILD NEW_OUTPUT
--export-archives` explicitly computes its two coarse two-k-point examples and saves
archives on the subgroup root (one or two ranks in this demo). Four SCF iterations are an ownership demo, not physical
convergence; this optional sample was not run in the present no-new-GS campaign.
Export copies full native orbital buffers including halos and creates reconstructed
density/chunk buffers; orbital_block controls disk chunks, not peak memory.

Regular k coordinates follow SALMON `lattice.f90:init_kvector`: `(i+0.5)/n-0.5`
with zero-based i, x fastest, then y, then z, and uniform weight1/nk. This is the
unshifted non-Gamma-centered full mesh generated by immutable `DFT.input_text`,
not a reconstruction of arbitrary user k tables or symmetry-reduced grids.
Density sums `weight[k]*occupation[k,b]*abs(u[k,b])**2` in native k/band order.
Synthetic2/3/6/8-point tests preserve each k/band, coordinates, weights and electron
count; Gamma output files remain byte-identical to commit b66b97fd. Custom/reduced-k retained owners and native GS-handle restoration remain future.
Native GS/2 import retains its own converter/runtime/cell restrictions: being
exportable does not guarantee compatibility with every native importer.

### MPI export ownership and failure behavior

All subgroup ranks enter matching `export_archive` calls. k/orbital/spatial and
combined partitions use explicit native lower bounds and physical density bounds,
not an assumed rank numbering. Replicated occupations/eigenvalues must match;
physical orbital tiles must cover global k/band/grid exactly without overlap.
The subgroup root writes GS/2 and returns GroundStateArchive; workers return None.
Only the root needs a writable destination. Explicitly open the archive on other
ranks later if they share that filesystem. Calling `gs.close()` remains collective.

Root receives one tile per collective to avoid retaining every payload at once;
it still holds the full global complex wavefunction, local copies, reconstructed
density and writer buffers. There is no arbitrary memory cutoff. A real allocation
or write failure is broadcast to all ranks without mutating native GS owners;
new incomplete archive directories are cleaned, existing outputs preserved.
Changing MPI topology during native restoration is not implemented by this export.
General failed MPI communications/process recovery is outside this contract.

Synthetic actual-MPI4 tests cover k/orbital/spatial/mixed partitions (including
uneven band splits), serial-file parity, failed bounds/replicas/coverage/metadata/
options/allocation/write operations, retained source arrays and successful retry.
This does not validate actual SALMON GS export or SCF/RT physics. The runner is
`tools/test_retained_gs_archive_mpi.py NEW_OUTPUT` under mpiexec -n4, also installed
under salmonpy/_runtime_source/tools; it never loads a native GS/runtime.
