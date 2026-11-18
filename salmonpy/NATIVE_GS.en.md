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
state with a host MPI subgroup. Both Gamma and multiple k points are supported. Raw legacy-file imported owners lack a complete model
and cannot use this entry yet; GS/2-imported owners use the context described below; complete native SCF-state restoration remains future work.
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
count; Gamma output files remain byte-identical to commit b66b97fd. Custom/reduced-k retained export and complete native SCF-state restoration remain future.
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

### Referenced input assets on existing-GS import

`session.load` retains/copies input-referenced pseudopotentials, coordinate files
and `file_kw` k-point tables with their relative subdirectories. These bytes are
verified on every subgroup rank and retained for in-memory GS→TDDFT models.
Quoted `!` in filenames is preserved while actual input comments are ignored.
Missing referenced assets fail admission before native initialization; accepted
GS/TDDFT owners remain unchanged. An unreadable required file is not a numerical
convergence warning: native initialization cannot consume it.

The shared existing-GS comparison runner now accepts `--external-k`:

```sh
mpiexec -n 4 python tools/test_memory_gs_mpi.py BUILD EXISTING_GS NEW_WORK --size 2 --split orbital --external-k
```

It stages a scratch copy of the existing input/GS, replaces a regular unshifted
source mesh with the identical explicit table, and compares against the original
file route. It does not run SCF or create a reference GS. Source inputs are unchanged;
only diagnostic summaries should be collected from its scratch output. The optional
runner is restricted to regular meshes without explicit shifts/Gamma centering;
this is not a restriction added to `session.load` itself.
PZ2-k, Libxc-PBE8-k and TB-mBJ8-k checks cover internal MPI2, independent/sibling
states, direct no-replay restart and missing-asset rejection while RT is active.
The runner is also bundled under salmonpy/_runtime_source/tools.
Imported `calculation` remains None; raw legacy model binding, full GS/SCF-state
restoration and physical source provenance remain separate. GS/2-imported owners
can re-export using the retained source context described below. File integrity and trajectory agreement do not certify physical
SCF convergence.

## Captured declared-input metadata

Install the optional `metadata` extra (f90nml) on the management machine before
moving to offline compute nodes. NumPy remains the only required dependency;
ordinary import/run/load does not import f90nml or acquire anything.

```python
# Inside an existing session; input_metadata itself is local, not collective.
metadata = gs.input_metadata()
# All ranks still call collective close as usual.
gs.close()
record = metadata.as_dict()  # detached JSON-compatible data, no native handle
cell = metadata.settings["cell_bohr"]
```

`GSInputMetadata` is immutable; `as_dict` returns a detached editable copy.
It retains user material_id, captured input/asset SHA256 values and a source_id.
source_id identifies those declared bytes, not physical material equivalence,
binary authenticity or checkpoint compatibility. Reformatting input or changing
its file names can change it. The native identity and user material ID are unchanged.
`gs.input_metadata()` can be called after close: captured text/assets are retained,
without reading source files, native calls, MPI or implicit recomputation.

The description normalizes explicitly declared au/A_eV_fs cell lengths to bohr,
ordered species and supported inline/external Cartesian/reduced atomic records
to fractional positions. It reports explicit grid counts/spacing, regular k-mesh
settings or external table coordinates/input weights, functional settings and
species-indexed potential hashes/local channels. Native-normalized k weights and
omitted defaults are not inferred; ambiguous/unsupported geometry remains None
with notes. Unknown source units do not silently become atomic units.
This is a declared-input description, not a query of the effective native model.

Imported evidence retains unknown convergence/provenance; computed owners retain
the convergence flag already obtained during run. Reading metadata does not run
SCF, recertify a GS, enable uncomputed forces, set imported calculation, authorize
native restoration or enable imported-owner portable export. No new restrictions
are imposed on ordinary load/run. Unsupported optional descriptions do not modify
GS state; callers may catch dependency/input errors independently of computation.

For custom hosts, `salmonpy.gs_input.describe_gs_input(text, assets=relative_name_to_bytes,
material_id="film")` uses the same description without a runtime or filesystem.
`python -m salmonpy.examples.memory_gs BUILD EXISTING_GS NEW_WORK --describe-input`
prints the root's snapshot after GS close. The MPI example agrees on root errors
before entering further collective operations. The reusable existing-GS test runner
supports `--describe-input` too. Optional analysis requires f90nml on the rank(s)
that call it; import or help does not initialize MPI.

## Initialize a retained GS from portable arrays

`session.load_archive` explicitly connects GS/2 to the existing native loader.
Prepare a probed GammaRestartConverter or SerialRestartConverter and its paired
serial SalmonExecutable beforehand. No implicit compiler, download or SCF runs.
All subgroup ranks call with matching selections; only root reads/converts the
archive, and the new converted directory must then be visible to subgroup ranks.

```python
with backend.open(context) as session:
    gs = session.load_archive(
        "portable-gs", "new-converted-bundle",
        converter=converter, executable=paired_serial_engine,
        pseudopotentials=potentials, material_id="silicon")
    origin = gs.archive_origin  # immutable hashes/lineage; valid after close
    model = gs.tddft_model(dt=0.003, max_steps=100)
    gs.close()
    state = model.new_state("response")
    j = state.step([0.0001, 0, 0])  # transverse A/c -> charge current
    state.close()
    model.close()
```

Omitting material_id inherits retained-export provenance's material_id when present,
otherwise "material". An explicit ID is honored; archive_origin retains the source
archive ID and convergence claim separately from current import convergence(None).
Manifest/converter/paired-engine/shared-runtime hashes bind the route. Input-derived
metadata remains available through input_metadata; calculation remains None and
forces uncomputed. Source archives and earlier owners are not overwritten.

This initializes orbitals/occupations in a new native owner. Archive density,
eigenvalues and energy are not serialized native restart fields; initialization
reconstructs native potentials/density. Mixing history, native solver caches and
in-progress GS iterations are not restored. Conversion errors are agreed before
native calls; new incomplete conversions are cleaned. A completed converted bundle
is retained for the caller, including if later native admission fails. Existing
bundle paths are preserved. Native STOP/process/MPI failures retain existing limits.

Converter restrictions remain: periodic unpolarized fixed ions, closed-shell even
electrons, grid>=8, PZ/Libxc-PZ/PBE with a supporting shared runtime; Gamma converter
is exact-Gamma, multi-k SerialRestartConverter uses aligned positive orthogonal
cells. Conversion is root/serial, while native initialization supports compatible
k/orbital/spatial subgroups. This is not arbitrary RT/SCF checkpoint repartitioning.
The paired serial executable is a format identity here, not the running solver.

Recorded incompatible Fortran record-marker/integer widths are rejected before
conversion. Missing compiler/source attestation or other flag differences warn and
proceed with the explicitly chosen native loader. Legacy receipts often lack a
compiler hash; no artificial attestation is invented. General native I/O ABI safety
across unvalidated compilers/platforms is not promised. Native record environment
must agree across ranks; converted receipts/contents are verified before loading.

The sample `salmonpy.examples.portable_gs.current_from_archive(session, ...,
A_steps=..., dt=...)` returns j for A/c at dt,2dt,... and closes its own owners.
Keep a single session for repeated calls; each conversion directory is new and
remains caller-owned. Prepare converter tools on management nodes and stage their
validated deployment records/binaries before using offline compute nodes.

Gamma array imports also use an equivalent explicit Gamma table so native DFT
allocates complex orbitals; nonzero imaginary components are preserved rather
than entering the default real-orbital Gamma branch. The converted k-table hash
is included in its receipt, with this route marked in archive_origin. Synthetic
complex Gamma tests pass under orbital/spatial MPI, including exact readback and
independent/file/current/direct-resume controls.

## Re-export a GS/2-imported owner

A live owner returned by `load_archive` can use the same collective
`export_archive` entry. Raw legacy `load` owners still need a complete physical
model and are not silently promoted from partial input descriptions.

```python
with backend.open(context) as session:
    gs = session.load_archive("source-gs", "new-native-bundle",
        converter=converter, executable=paired_serial_engine,
        pseudopotentials=potentials, material_id="silicon")
    saved = gs.export_archive("new-gs-array-archive")  # root returns archive, others None
    # GS and any existing RT models remain unchanged.
    gs.close()
```

Only small source metadata/eigenvalue/k/geometry arrays and immutable potentials
are retained during load. Source wavefunctions or a live source archive object
are not cached. Re-export can work after source-archive removal, while the native
GS owner remains live. It reads native orbitals/occupations and removes halos.
Density is reconstructed with source sampling weights. Geometry/grid/XC/potential
references, k coordinates/input weights, reference eigenvalues and energy are
carried from the source archive rather than claimed as freshly computed results.

The optional provenance variant `salmonpy.retained-gs-reexport/1` explicitly labels
energy_source/eigenvalue_source/sampling_source as source_archive, records the
source manifest SHA and inherited convergence claim, and leaves current converged
None. Source/reference values do not imply SCF convergence or native cache
restoration. This provenance variant requires a reader supporting it; existing
GS/2/computed-export provenance is unchanged. calculation stays None and forces
remain unavailable. Complex reads copy real/imaginary components without
arithmetic, preserving signed-zero bits and rejecting nonfinite values collectively.

Nonuniform weights and irregular k sequences are preserved rather than inferred
as a regular mesh. Existing computed-DFT export keeps its arithmetic order. The
sample current_from_archive(..., reexport_to="new-array-archive") optionally saves
before creating RT states. Repeated save/load initializes new models, not a transfer
of an existing RT checkpoint to another model/topology. Full SCF mixing/cache/
iteration continuation remains separate.
