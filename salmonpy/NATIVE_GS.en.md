# Retained native ground states

> Current build layout (2026-11-19): ordinary XC and fast PR HSE both build
> from runtime/sources/solver. Below, standard/ordinary identifies the API
> capability profile, not a second maintained source tree. Default builds disable HSE.

[日本語](NATIVE_GS.md) · Topic: `retained_native_gs`

With a standard runtime, salmonpy.native_gs keeps several DFT ground states in one process.
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

Python scope is periodic full-k unpolarized fixed-ion PZ/Libxc-PZ/PBE and HSE06
with a declared fast-PR compute runtime. HSE uses k-only MPI by default; explicit Gamma orbital GS is described below. SOC/+U/DC GS,
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

### Collective array requests

Call `gs.array(quantity)` on every rank in the GS subgroup, selecting the same
supported quantity string. Invalid types, unknown quantities, unavailable imported
forces and differing rank selections raise `InputError` collectively before native
array access. A corrected request can be retried; retained GS and active TDDFT
owners are unchanged. This cannot recover a missing rank or mismatched API call
order. Numerical values and supported quantities are unchanged.

### SCF continuation groundwork

New standard-runtime builds retain the SCF convergence work arrays (previous
density and local potential), electron normalization, no-mixing counter and final
iteration/residual/convergence/ion-update values inside each owner's mixing state.
The existing entry still resets this progress for each fresh SCF call and keeps
the original arithmetic and convergence policy. Existing GS files, GS/2 import
and RT checkpoint formats are unchanged. This is internal groundwork: no public
SCF resume API or complete SCF checkpoint is available yet. The two previous-field
arrays now remain live with the GS owner until reset/close. HSE source profiles
are unchanged by this standard-only slice.

## Fast-profile existing GS import

The standard SCF example above remains supported by the standard runtime. A fast
build declaring `gs_abi` protocol1 also supports existing-GS import, without SCF:

```python
context = SimpleNamespace(solver_comm=MPI.COMM_WORLD, closed=False, _dependents=0)
backend = NativeGSBackend(fast_build_record, shared_work_root,
                         decomposition=(MPI.COMM_WORLD.size, 1, 1, 1, 1))
with backend.open(context) as session:
    a = session.load(bundle_a, material_id="film-a")
    b = session.load(bundle_b, material_id="film-b")
    density = a.array("density")  # independent read-only copy
    occupations = b.array("occupations")
    a.close()                     # b remains readable
```

The runtime record's declaration is required; older fast libraries are not
silently accepted or rebuilt. Inputs describe compatible existing wavefunctions;
full/reduced k-point initialization keeps each model's symmetry scratch separate
and copies the required sym.dat. Import disables unused RT factor-learning and
pre-SCF controls, preserving XC and HSE screening. All calls remain subgroup
collective on the main thread; one attachment lifetime/process as above.

Fast imports expose density, occupations and orbitals. Their info energy_hartree,
residual, tolerance and converged are None, not physical zero. Eigenvalues and
forces were not computed. For import-only records session.run remains unavailable; imported owners cannot export computed results. tddft_model
requires the additional memory-RT capability described below; the standard profile retains its existing features. No convergence
claim follows from importing a synthetic or unlabelled GS.

The existing sample supports this path with no new sample file:
`mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD NEW_OUTPUT --import-bundle GS_A --second-bundle GS_B`.
Each state keeps its material_id; output JSON records unknowns as null. Testing
uses synthetic inputs and validates state independence, not physical convergence.


### Retained GS to fast-PR HSE response

A newly built fast runtime with `gs_abi.memory_rt == 1` and
`sp_fast_hse_memory_gs_protocol() == 1` connects `gs.tddft_model(dt=..., max_steps=...)`
to the fast-PR HSE registry. This slice requires HSE06, atomic units, a compatible
periodic scalar GS and k-only one/two-rank groups. Existing standard GS→TDDFT
support is unchanged; ordinary PZ/TB-mBJ/Libxc RT is available in the HSE-disabled fast profile below.

```python
with backend.open(context) as session:
    gs = session.load(bundle_a, material_id="film-a")
    model = gs.tddft_model(dt=0.001, max_steps=32)
    gs.close()  # The RT model owns a detached value copy of initial orbitals.
    a = model.new_state("point-a")
    b = model.new_state("point-b")
    try:
        j_a = a.step([0.0, 0.0, 1.e-5])
        j_b = b.step([0.0, 0.0, -1.e-5])
    finally:
        a.close()
        b.close()
        model.close()
```

No SCF, GS restart copy/output, or accepted-step replay is used for this handoff.
Static potentials/symmetry and RT resources are initialized normally. Each model
retains its own initial orbitals/occupations, so this costs memory; it is not a
zero-copy interface. Later source GS close/removal and sequential GS imports do
not replace existing RT model settings. Keep the session/host communicator alive
until every RT state closes. Model close prevents creation while live children
remain usable. Calls are collective, main-thread and sequential.

A is Cartesian A/c; returned j is charge-current density. Full-k inputs allow
x/y/z and rotating drives. Reduced symmetry must preserve the actual drive; the
example uses z for the tested C2z fixture, not a universal polarization limit.
Direct rank-local checkpoints use the existing same-runtime/topology contract.
The entry initializes from GS arrays, not an arbitrary RT/SCF checkpoint.

The existing sample adds `--current-steps 16` to its `--import-bundle` mode.
`tools/test_memory_gs_mpi.py BUILD FIXTURE NEW_OUTPUT --size 2 --fast-hse`
uses four ranks to compare memory and existing-file initialization. Synthetic
Taylor/learned PT-CN full/reduced controls test currents/energy, sibling ownership,
GS close, new states after source removal, interleaved GS import and direct
restart without replay. These are implementation checks, not physical convergence.


### Ordinary XC on the fast source profile (HSE disabled)

The same prepared fast-PR source can now build a normal CPU host-MPI library:

```sh
python tools/prepare_hse_rt_state.py NEW_SOURCE --explicit-routing
python tools/build_embedded_tddft.py NEW_SOURCE NEW_BUILD --jobs 4 \
  --cmake-arg=-DUSE_HSE=OFF --cmake-arg=-DUSE_LIBXC=ON
```

Supply the MPI/compiler/BLAS/Libxc locations for your environment explicitly.
This profile does not require HSE/FFTW dependencies. Its build record keeps
`salmonpy.embedded-tddft/1`, declares `source_profile="fast-pr"`, `ordinary_rt_abi=1`
and `gs_abi/1`. New HSE-disabled builds also declare the run mode. Use the existing NativeGSBackend session.load and
gs.tddft_model APIs or NativeTDDFTBackend file/registry APIs. PZ, TB-mBJ, Libxc-PZ
and Libxc-PBE, periodic full-k scalar fixed-ion middlepoint A/c response are
supported. Material IDs persist; returned current is electric charge current j.

HSE-disabled builds expose session.run for SCF, energy, forces and eigenvalues.
Nonconvergence warns and retains the last iterate. Imported metadata remains None.
Computed GS owners connect to explicit export_archive; this validation inspected
assembled export arrays in RAM without writing a GS archive. Source GS close, new independent states and direct
RT restart do not require GS-file copying or step replay. Model seed copies
consume memory; keep the session and communicator alive for all RT children.
Sequential subsequent GS imports suspend/resume RT inputs/timers. Main-thread,
collective and one attachment-lifetime constraints remain.

The ordinary registry reuses the standard embedded-TDDFT source, with compile-time
selection of the fast GS seed provider. Fast stencil and XC scratch already belong
to numerical states, so legacy global-workspace switching is omitted there.
For file initialization, parser snapshots capture API dt/capacity/MPI decomposition
after worker configuration, matching the standard reader's order. The existing
fast HSE ABI and standard runtime behavior are preserved. This does not yet
unify ordinary XC and HSE registries in one GS session or remove the old source
tree. New fast-HSE compute records are described below; final tree consolidation remains.

The retained_gs sample's import mode also accepts these records with
--current-steps 16. Native k/orbital/spatial MPI2 response, direct restart and
different-XC GS interleaving are tested with synthetic fixtures. Old-runtime
comparisons distinguish exact PZ/TB-mBJ agreement from tiny Libxc differences
(max5.5e-20au in the measured controls). Resume/fresh-state comparisons remain
exact. This is execution/numerical parity, not physical convergence.

Computed PZ/Libxc-PZ/PBE controls used orbital and spatial MPI2: independent
owners matched exactly. PZ/PBE energy and text reports matched the previous fast
executable. Four-iteration unconverged PBE still differs from the old standard
profile by about 3.4e-7 Hartree; no standard-profile identity or physical
convergence claim is made.

Computed GS owners use the same tddft_model entry. With an HSE-disabled,
Libxc-enabled build record, this example performs four SCF iterations, closes
the first GS owner and advances its independent RT states for 16 steps:

```sh
mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD_JSON NEW_OUTPUT --current-steps 16
```

This tiny example is not physically converged. SCF nonconvergence only warns.
It writes GS scalar summaries and charge currents to JSON, without GS wavefunction
or restart output. Material IDs persist; two independent states receive an A/c
drive with x/y/z components. Existing --import-bundle usage remains available.

Integration controls retain the RT model while another GS is computed, using
PZ/Libxc-PZ and Libxc-PBE with orbital/spatial MPI2. Currents after source GS close
match independent execution exactly. A direct checkpoint restored halfway through
16 steps continues exactly, without replaying time propagation.

### Computed GS with fast-PR HSE

New fast-HSE builds declare run in gs_abi.modes. Older import-only builds are not
silently upgraded. Pass DFT(..., xc="hse06") to NativeGSBackend session.run with
decomposition=(size,1,1,1,1). Scope is periodic full-k unpolarized fixed-ion.
The existing GS→TDDFT model owns its seed after the source GS closes.

```sh
mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD_JSON NEW_OUTPUT --xc hse06 --current-steps 16
```

New HSE inputs explicitly set exx_pre_scf_threshold=0, so the requested SCF
iterations start in HSE rather than a PBE warmup phase. Nonconvergence warns and
returns the last iterate. Batch DFT.run HSE reports remain unconnected and direct
users to this native entry. No conventional-HSE fallback is used. Physical force
accuracy and convergence remain separate validation tasks.

Darwin HSE libraries reuse the ordinary TDDFT linker policy to bind the selected
BLAS, avoiding NumPy's incompatible complex-return ABI. Explicitly build a new
runtime; fixed previous libraries are not modified.

### Batch energy/forces without restart arrays

Batch DFT.run keeps its existing save-by-default behavior. Set save_restart=False
per call to obtain energy and forces without GS restart arrays:

```python
gs = calculation.run(engine, "new-output", scf=SCFSettings(max_iterations=4),
                     save_restart=False)
print(gs.energy_hartree, gs.forces_hartree_per_bohr, gs.restart_saved)
```

Batch scope is PZ/Libxc-PZ/PBE. Four iterations test execution, not physical
convergence. Unconverged results warn and remain usable. restart_saved reports
whether all three restart files actually exist. Energy/force access works without
them; band/TDDFT initialization, warm start and GS archive export require a result
with restart arrays. Legacy manifests without this flag retain saved semantics.
If a runtime ignores output suppression, warn and record actual output rather
than turning a completed calculation into a failure.

The fast source now emits gs_status.dat (SALMONPY_GS_V1, Hartree) and gs_forces.dat
(SALMONPY_FORCES_V1, Hartree/Bohr), preserving SCF/force operation order. Native
HSE and CLI produce these reports too; batch DFT.run HSE selection still requires
fast-PR provenance integration, so use NativeGSBackend for HSE. The existing
dft_si sample accepts --no-restart and --max-iterations for this workflow.

### NLCC migration before solver-source consolidation

The fast source now carries the standard profile's NLCC Libxc energy weighting
and core-force correction. The valence-plus-core density supplied to Libxc also
weights its energy density. Core forces differentiate the actually sampled core
table using each material's species, boundary configuration and XC potential.
Fast HSE projector derivatives, owned state and existing spin handling are retained.

V1 cutoff behavior and V2 full support remain distinct. V2 uses the same table
interpolation and periodic-image extent for orthogonal/skew cells and long cores;
its cancellation-sensitive force integral uses compensated summation. IEEE
inquiries remain scalar and only used V2 table entries are inspected, not padding.

Prescribed-density energy/force variation and OMP/MPI partition checks are
implementation tests, not material convergence certification. They do not attach
validated NLCC-force capability records to an ASE Calculator runtime. Existing
capability gates remain. No-NLCC workflows retain regression checks. Removing
legacy sources and changing default preparation still require the remaining
workflow/provenance validation; this migration alone does not complete unification.

## Earlier ownership stage (superseded by the active SCF section below)

The shared fast-source mixing state now owns previous density/local potential,
electron normalization, convergence flag/count, iteration/residual and ion-update
state. Assignment deep-copies its numerical arrays. Arithmetic and update order
are unchanged. Existing fresh-SCF calls reset this progress. This establishes
ownership for checkpointing, not a completed public no-reinitialization SCF
continuation or complete portable checkpoint. Mixing/CG/potential/cache continuation
remains next. Execution order is SCF persistence→shared ownership/recovery/native
ABI→NLCC/units. Fast-HSE norm diagnostics and a physical grid campaign are deferred
at the user's request. No GS datasets are collected by the validation probes.

## Active SCF pause and in-memory snapshots

New explicitly built CPU runtimes expose SCF protocol1 for periodic,
unpolarized fixed-ion PZ/Libxc-PZ/PBE (with or without explicit symmetry), and fast-PR HSE06
(full k or the reduced-k scope below). Existing run() remains a
full SCF call. A steps budget pauses only at an accepted iteration boundary:

```python
from salmonpy.dft import SCFSettings
with backend.open(context) as session:
    with session.run(calculation, scf=SCFSettings(max_iterations=100),
                     material_id="film-A", steps=5) as state:
        with state.snapshot() as saved:
            state.continue_scf(steps=10)
            state.restore(saved)  # zero SCF iterations, no initialization/replay
            state.continue_scf(steps=20)
        print(state.scf_status)
```

The configured total cap/tolerance/model stay fixed; steps is a per-call budget.
Convergence may finish earlier; reaching the total cap without convergence warns
and retains results. A deliberate pause is not reported as failed SCF. Forces
are available after the SCF segment finishes, not at an intermediate pause.

NativeSCFSnapshot retains independent numerical copies in Fortran memory:
orbitals/CG, occupations, mixing history/rates/progress, density/spin/potentials,
Poisson/reciprocal/XC workspace and effective post-initialization inputs. Original
declared model inputs remain separate for GS→TDDFT seeds. Continue bypasses
initialization, initial no-mix work and convergence reset. Restoring performs no
Hamiltonian solve or SCF replay. scf_status reports initializations and executed
main iterations; explicit branch calculations after restore add their own work.

All calls are collective, in identical order on the internal MPI group, on its
main thread. A snapshot belongs to the same live GS owner/build/decomposition.
Close snapshots explicitly while the owner is live; closing the owner frees all
remaining numerical copies. Complete copies consume memory. This entry is not a
disk checkpoint, pickle format or process-restart/repartitioning mechanism. Text
logs/output positions are not rewound. Invalid/mixed snapshot selections reject
before numerical writes; general OOM/initialization/process-failure atomicity is
not established. SOC/+U/DC/ionic-motion continuation is not connected yet. Fast-PR HSE06
requires its compute-capable HSE runtime and k-only MPI decomposition.
Old runtimes remain usable for ordinary run/load; new APIs require a matching
explicitly built runtime and never trigger an implicit build.

```sh
mpiexec -n 2 python -m salmonpy.examples.scf_resume BUILD_JSON NEW_OUTPUT
```

The sample performs a small execution check with GS restart output disabled.
It compares continued/restored numerical values in RAM and saves scalar evidence
only; it does not collect GS wavefunction/density datasets.

### Fast-PR HSE06 SCF snapshots

The same run(steps=...), snapshot/restore and continue_scf entry supports full-k
CPU fast-PR HSE06, including Gamma and internal k-MPI. Select
`decomposition=(comm.size, 1, 1, 1, 1)` and a matching HSE build; the ordinary
HSE-disabled runtime cannot execute this profile. There is no slow HSE fallback.

Snapshots also retain the explicit owner's ACE/damped operator, cached
source/action/occupation, factor/PT-CN history arrays and flags, exchange energy,
support decisions, actual mixing age and PBE pre-SCF phase/count/first iteration.
The owner's FFT/backend kernel stays live; it is moved aside during numerical
copying, never intrinsically copied, recreated or evaluated by restore. System
owner pointers are rebound to the same live context, not to the saved payload.
A subsequent natural PBE→HSE transition still builds the target exchange as in a
continuous calculation; that is future work, not restore-time replay.

A target-HSE iteration cap without convergence warns and retains the last
iterate, including after PBE pre-SCF; it is not labelled converged. An unfinished
PBE-only pre-stage is still not a target-HSE result. SCF thresholds in the small
phase-transition tests are control-flow settings, not physical certification.
The snapshots remain memory-only, same owner/build/decomposition. Reduced-k
HSE is covered below; SOC/+U/DC, ionic motion and
process/disk restart remain outside the current entry.

```sh
mpiexec -n 2 python -m salmonpy.examples.scf_resume HSE_BUILD_JSON NEW_OUTPUT --xc hse06
```

Native controls compare density, occupations, orbitals, eigenvalues, forces and
metadata to both continuous new and fixed preceding runtimes. Cuts before/after
the PBE→HSE switch retain mixing history; a separate tight-threshold test confirms
nonconvergence warns. Interleaved different-material snapshots remain independent.
GS→16-step TDDFT currents also match exactly. These are execution regressions;
they do not establish physical convergence. No GS datasets are collected.

### Common checkpoint binding

`state.checkpoint_capability` declares native-scf-memory/1 for this live owner.
`state.checkpoint_target()` binds save; `state.checkpoint_target(snapshot)` also
binds direct restore. Use execute_checkpoint from CHECKPOINTING.en.md; the scf_resume
example now exercises this common entry. Payloads stay local Fortran memory.
The binding does not execute a save or SCF iteration. New matching runtimes add
independent numerical candidate staging; after restore_checkpoint's ready vote,
publish each returned candidate before continuing (CHECKPOINTING.en.md).
Old runtimes reject protected staging, with no silent in-place fallback.
Multi-component direct callbacks remain non-atomic.

### Explicit symmetry and reduced-k HSE SCF

`DFT(..., symmetry=SymmetryOperations(operations, directions="yyn"))` attaches
explicit fractional rotations and translations. operations has shape (n,3,4):
the first three columns are the rotation and the fourth is the fractional shift.
Data are copied and immutable. No symmetry discovery, cell standardization or
network lookup occurs. Native fast-PR admission retains the existing fractional
signed-permutation subset and checks that cell/grid/atoms/k stars are compatible.
The user chooses a subgroup preserving the actual crystal and, for later RT,
the applied field; this does not impose a global z-polarization restriction.

```python
from dataclasses import replace
import numpy as np
from salmonpy.dft import SymmetryOperations
ops = np.zeros((2, 3, 4))
ops[:, :, :3] = np.eye(3)
ops[1, :2, :2] = -np.eye(2)  # C2 about z; identity is also retained
calculation = replace(calculation, symmetry=SymmetryOperations(ops, directions="yyn"))
```

`input_assets` includes generated sym.dat and pseudopotentials. Native/batch
inputs write these same owned bytes; MPI checks asset content before entry.
Symmetry hashes enter model/warm-start identity and assets follow GS→HSE TDDFT.
Default symmetry=None preserves existing input text and calculation paths.
Reduced-k native SCF currently uses CPU fast-PR HSE, full q, one/two internal
k ranks. Pause/continue/snapshot, direct restore and protected numerical staging
use the same APIs as full k. Explicit symmetric PZ/Libxc-PZ/Libxc-PBE SCF uses the same APIs.
Internal k-MPI2 and spatial-MPI2 controls verify interleaved-material isolation,
protected restore and exact uninterrupted parity.
Symmetry/density maps are owned by each GS context and swapped around continuation;
immutable initialization tables stay with their live model, not a new solver clone.

```sh
mpiexec -n 2 python -m salmonpy.examples.scf_resume HSE_BUILD_JSON NEW_OUTPUT --xc hse06 --symmetry c2z --protected
```

C2 reduced SCF/restore orbitals, occupations, energy/density and 16-step z-driven
RT match fixed preceding continuous results exactly; interleaved materials use
different cell/grid settings. The intentional symmetric-force units change is
compared after reconstructing the old basis, as explained below. These are execution tests, not physical convergence certificates.
A fixed-grid diamond-Si check separately SCF-converges full8/reduced6 k points:
energy agrees within1e-6 Hartree and density within rtol1e-4/atol1e-7 bohr^-3.
The comparison exposed an existing units defect: force averaging converted
Cartesian force into direct-lattice vector components but failed to convert the
average back. Identity/C2/skew-cell probes reproduce the failure on the old
routine. The maintained Fortran now returns `Amat * average(Ainv * F)`.
**Symmetric force values intentionally change to Cartesian Hartree/bohr.**
Nonsymmetric forces, energies, densities and fixed-ion currents are unchanged.
Old symmetric force output is compared only after reconstructing its Cartesian
basis in the regression test; that conversion is not an application fallback.

The same fixed-grid, SCF-converged full8/reduced6 Si comparison now also passes
the1e-6 Hartree/bohr force criterion (largest difference about1.15e-7).
This checks symmetry equivalence and units, not complete physical force, grid/k,
finite-displacement or geometry-optimization certification. NLCC/ASE capability
gates are not relaxed. Historical outputs are preserved, not rewritten.

A new runtime query sp_gs_cartesian_force_protocol()==1 declares this basis.
For older symmetric runtimes without it, Python warns, keeps values available,
and labels their basis unverified rather than falsely claiming Hartree/bohr.
Use a matching corrected runtime for Cartesian symmetric forces. There is no
implicit runtime build or Python numerical conversion. Reduced-k GS archive
export, disk/process restart and topology changes remain separate work.

The small developer probe is `python tools/test_symmetry_force_units.py NEW_OUTPUT`
(requires an explicit Fortran compiler; it compiles/audits the actual module, not
a new GS). It is not a prerequisite for ordinary prebuilt-runtime use.

### Disk SCF checkpoints and fresh-process continuation

A matching new runtime supports `state.save_scf(new_directory)` and loading in a
fresh host with `session.run(calculation, scf=scf, material_id=..., restart=directory)`.
This differs from reusable converged-GS `export_archive()` output. The runtime
SHA256, material ID, input/pseudopotential/symmetry assets, rank count and actual
MPI decomposition must match. This is not a changed-input warm start or repartitioned restart.

```python
# Saving process
with session.run(calculation, scf=scf, material_id="Si", steps=3) as state:
    state.save_scf(checkpoint_directory)  # New directory; never overwrite.
# Fresh process, same calculation/scf and decomposition
with session.run(calculation, scf=scf, material_id="Si", restart=checkpoint_directory) as state:
    assert state.scf_status["executed_iterations"] == 0
    state.continue_scf(steps=5)
```

Fortran writes rank-local stream arrays, density-mixing history/rates, accepted
iteration/no-mix progress, XC phase, and fast-HSE ACE/history/occupation/cache.
Python authenticates versioned metadata and SHA256; the manifest is written only
after every rank succeeds. Incomplete output is not a valid checkpoint. Use the
calculation-results area for saved GS data.

Load validates all rank files and model inputs before creating a context. Static
communication/grid/pseudopotential/Ewald and FFT/Libxc resources are prepared;
arrays are read into detached candidates and published only after collective
success. Initial guessed orbitals/density, mixing initialization, initial energy
evaluation and all SCF iterations are skipped. `executed_iterations` counts work
actually performed in the current process; prior counters remain in the manifest.
FFT handles, pointers and MPI communicator IDs are not serialized. The HSE FFT
kernel is prepared when the next natural iteration needs it, without replay.

PZ/Libxc-PZ/PBE and fast-PR HSE controls cover k-MPI2, ordinary-XC spatial-MPI2,
C2 reduced k and the PBE-prestage→HSE boundary: save at3 plus5 iterations in a fresh
process exactly matches continuous8. A corrupt rank file is rejected before
state creation; native read failure preserves the existing state. These are
execution regressions, not physical convergence certification. Nonconvergence
still warns and retains the result.

The existing scf_resume sample accepts `--save-scf NEW_CHECKPOINT` to save the
paused boundary. A separate invocation uses `--restart CHECKPOINT`; keep material
settings, `--xc`, symmetry, `--iterations` and MPI decomposition identical. The
accepted SCF position is restored from the file. Format: salmonpy.native-scf-disk/1;
optional ABI: sp_gs_scf_disk_protocol()==1. Different binaries, changed topology
and portable GS conversion remain outside this format.

### Residual history and unconverged-SCF plots

`save_scf()` also writes diagnostics/residuals.csv, diagnostics/residuals.png and
diagnostics/report.json. The logarithmic figure shows residual detail and a second panel overlaying
the residual and threshold, and explicitly labels unconverged results. CSV columns are accepted iteration,
residual, threshold and PBE-prestage flag. Initial no-mixing iterations do not
evaluate this residual and are excluded. The residual follows the chosen native
SCF convergence metric; it is not a physical convergence certificate.

`state.array("scf_history")` returns the same read-only4×iteration-count array.
Fortran retains the history across memory snapshots and disk restart, so a report
can show both sides of a resumed run. Use `state.write_scf_report(new_directory)`
for a report without saving another checkpoint, including after continuation.

Matplotlib is optional. Missing plotting support or a plotting failure only
warns: the checkpoint and CSV remain available. A plotting failure must not be
reported as a failed SCF calculation.

### Common disk execution and protected restore

`state.scf_disk_capability` declares the format, actual MPI decomposition and
material/input/runtime identity, without file I/O. Bind
`state.checkpoint_target(directory=path)` to execute_checkpoint/restore_checkpoint.
Choose either the memory snapshot argument or directory, never both.

```python
from salmonpy.checkpointing import execute_checkpoint, restore_checkpoint, CheckpointState
# a,b are independent SCF states on the same solver_comm; matching rank call order.
groups = [{"name": "scf", "ranks": list(range(comm.size)), "point_ids": ["A", "B"]}]
targets = {"A": a.checkpoint_target(directory=checkpoint_a),
           "B": b.checkpoint_target(directory=checkpoint_b)}
saved = execute_checkpoint(comm, groups, states=targets)
# After advancing the original states, restore the saved boundaries.
live = CheckpointState(states={"A": a, "B": b})
prepared = restore_checkpoint(comm, groups, states=targets,
                              saved=saved["compatibility"], live=live)
for candidate in prepared["states"].values():
    candidate.publish()  # Publish every candidate before resuming computation.
```

Staging authenticates files/inputs and reads detached numerical candidates; the
common API votes across participants. Handled staging failure discards candidates
and leaves originals available for continuation. Recreate checkpoint targets for
a retry after failure. Formats and I/O remain backend-owned under the existing
common contract.

For one owner, `state.restore_scf(directory)` uses the same staged implementation.
Direct calls across multiple owners are not atomic as a whole: use
restore_checkpoint for coordinated protection. Save/restore performs no SCF
iterations or resource cloning. Fatal allocation/process/publication failure is
outside this rollback guarantee.

The existing scf_resume sample uses `--disk --save-scf NEW_CHECKPOINT --protected`
for the common disk path. Omitting --disk retains its existing memory behavior.

### Numerical-candidate cleanup

Prepared native SCF candidates support `candidate.close()` (equivalent to discard)
and context-manager use. Leaving the context before publish discards only the
candidate and preserves original numerical values. Closing an already published
or discarded candidate is a no-op; a closed candidate cannot be published.

Closing the Python GS owner first discards its Python-owned pending candidate.
Closing a memory snapshot first discards an unpublished candidate derived from
that snapshot. Host-code exceptions can therefore unwind these contexts without
leaving GS/MPI resources reserved. A snapshot does not discard unrelated candidates.

All MPI ranks must call close/publish/discard in matching order. No garbage
collector/destructor performs MPI operations. Reservations staged directly by an
external C host remain that host's responsibility; Python never silently discards
them, and native admission/refusal remains in effect.

### Changing k-MPI rank count at restart

A matching new runtime supports periodic scalar fixed-ion PZ/Libxc-PZ/Libxc-PBE
and fast-PR HSE06 SCF restart between k-only layouts; measured controls cover1→2 and2→1 ranks.
Default restart still requires the same decomposition. Opt in explicitly:

```python
# Same calculation/scf/material_id and runtime; open the new host MPI group with
# decomposition=(comm.size,1,1,1,1).
with session.run(calculation, scf=scf, material_id="Si",
                 restart=checkpoint_directory, repartition=True) as state:
    assert state.scf_status["executed_iterations"] == 0
    state.continue_scf(steps=5)
```

New manifests record each file's global k bounds. Only k parallelism may change:
keep cell/grid/k set/electron count/XC/mixing/tolerance/total cap/assets identical.
Admission validates format, runtime SHA256, material/input identity, contiguous
nonoverlapping old k ranges and all old rank-file hashes before creation.
Fortran reads overlapping donors and reindexes every wavefunction/CG k slice,
retaining density/mixing/Poisson/XC numerical state and residual history. The
continuation's parallel inputs bind to new live resources. Native staging also
checks complete target-k coverage before publication. Source files stay intact.

Immediately restored orbitals/density/occupation/eigenvalues/history exactly match
the saved boundary. Five subsequent iterations are compared with continuous8
from the same source-initialized state at rtol1e-9/atol1e-10. The default Gaussian seed includes the first local k index, so fresh
initialization with a different decomposition may produce different orbitals.
Initialization itself is unchanged; a different guess is not the same initial condition. These are execution controls, not physical
grid/k convergence certification.

Fast-PR HSE additionally reindexes dense ACE factors, source/action caches and
per-k factor/PT-CN histories in Fortran. Full/reduced-k and PBE pre-stage cuts
cover1↔2 ranks. Packed/distributed ACE and changed-runtime migration remain outside this entry;
HSE spatial restart uses the explicit Gamma entry below. Common checkpoint_target/protected restore still requires unchanged
MPI decomposition. Optional ABI: sp_gs_scf_k_repartition_protocol()==1 for ordinary XC, version2 for HSE; no implicit build.

The existing sample saves with `--k-parallel --save-scf NEW_CHECKPOINT` and loads
with a different rank count using `--restart CHECKPOINT --repartition`. Both
ordinary-XC invocations use the same2×2×2 k mesh; keep other inputs identical.

### Ordinary-XC orbital-only SCF repartition

The same explicit `run(restart=..., repartition=True)` entry also accepts pure
orbital layouts, with `decomposition=(1, comm.size, 1, 1, 1)` on both hosts.
The optional `sp_gs_scf_o_repartition_protocol()` version1 or2 ABI is required. New saves
record global orbital bounds; older checkpoints without these bounds retain
same-layout restore but cannot migrate orbital layout. All physical inputs,
assets and runtime must match. PZ/Libxc-PZ/PBE1↔2 controls compare exact restored
orbital slices, zero SCF bootstrap work and continued trajectories at the existing
rtol1e-9/atol1e-10. Fortran reindexes orbital/CG buffers and retains replicated
density, mixing, XC, Poisson and residual history. Real or complex orbitals use
the same native path. No Python numerical rearrangement occurs in the loader.

The sample accepts `--orbital-parallel --save-scf NEW_CHECKPOINT`, followed by
`--orbital-parallel --restart CHECKPOINT --repartition` on a fresh rank count.
Its orbital layout uses a fixed Gamma mesh unless explicit symmetry is requested.
Fast HSE orbital repartition uses the explicit protocol2 entry below. Changes between k-only and orbital-only layouts,
Ordinary spatial/mixed restart uses the entries below; changed-runtime migration remains unsupported.
Default and common protected restore still require identical MPI layouts.

### Fast-PR HSE Gamma orbital SCF: same-layout direct restart

A new explicitly built runtime declaring `sp_gs_hse_orbital_protocol()` version1 or2
accepts retained fast-PR HSE06 GS with `decomposition=(1, comm.size, 1, 1, 1)`
and Gamma-only sampling. The CPU full-support profile is periodic, unpolarized,
fixed-ion, orthogonal-cell, integer fully occupied spin pairs, with no symmetry,
localization or pair-screening approximation. These are the existing Gamma
exchange implementation's requirements. The k-only HSE path is unchanged.

`save_scf()` and `run(restart=...)` preserve the same runtime/input/MPI layout.
Wavefunction/CG, density/mixing, ACE and Gamma exchange workspace (including
its progress/gauge/source data) are owned by each state and serialized directly.
Backend procedure pointers and communicator IDs are not serialized; borrowed
communicators are rebound on restore. No SCF iteration or accepted ACE rebuild
is used to reconstruct the saved numerical state.

The existing sample uses `--xc hse06 --orbital-parallel --save-scf NEW_CHECKPOINT`,
then `--xc hse06 --orbital-parallel --restart CHECKPOINT` with the same rank count.
Use the protocol2 entry below when changing ranks. Multi-k orbital HSE,
Mixed layouts and orbital GS→RT handoff remain unsupported; HSE spatial GS
uses the explicit entry below.
Two-rank fresh-process continuation and a cut during the PBE pre-stage are
execution controls; they do not establish physical grid/k convergence.

### Changing fast-PR HSE orbital rank count

Use the explicit, immutable backend option `hse_orbital=True` on both saving and
restoring hosts, including one rank. It selects the same full-support Gamma
orbital exchange algorithm instead of interpreting a one-rank layout as the
default k-only algorithm. This requires HSE-orbital protocol2; migration also
requires orbital-repartition protocol2. Older default k/same-layout behavior
is unchanged. Every physical input and the runtime must remain identical.

```python
backend = NativeGSBackend(build_record, work_root,
    decomposition=(1, comm.size, 1, 1, 1), hse_orbital=True)
with backend.open(context) as session:
    with session.run(calculation, scf=scf, material_id="Si",
                     restart=checkpoint, repartition=True) as state:
        assert state.scf_status["executed_iterations"] == 0
        state.continue_scf(steps=5)
```

Fortran maps dense ACE factor columns and source/action/work/spatial-source/
previous columns by saved global orbital index. Replicated gauge, occupation,
phase and density-mixing histories stay intact. Communicators are rebound;
no accepted ACE evaluation or SCF replay reconstructs the checkpoint. Source
files remain unchanged. Fresh-process controls cover1↔2, PBE pre-stage1↔2 and
unequal two-atom/four-band2↔3 partitions. Restored orbitals are exact, bootstrap
SCF work is zero, and five new iterations match the same source-initialized
continuous8 at the existing rtol1e-9/atol1e-10.

The sample uses `--xc hse06 --orbital-parallel` on both sides; add `--save-scf`
on the source and `--restart CHECKPOINT --repartition` on the new rank count.
The constructor option is needed even when the source has one rank. Packed/
distributed ACE, tiled gauges, localized sparse support, multi-k orbital HSE,
Mixed layouts, changed-runtime migration and orbital GS→RT remain outside
this dense Gamma contract; HSE spatial GS uses the dedicated entry below. Common protected restore still keeps the same layout.

### Ordinary-XC spatial rank-count restart

A matching runtime declaring `sp_gs_scf_r_repartition_protocol()` version1 or2 accepts
periodic scalar fixed-ion PZ/Libxc-PZ/PBE SCF restart between pure spatial
layouts. Use `decomposition=(1, 1, nx, ny, nz)` with `nx*ny*nz == comm.size`
on each fresh host. Only the parallel block may differ; runtime, material ID,
physical inputs, assets, SCF cap and mixing settings must match.

```python
backend = NativeGSBackend(build_record, work_root,
    decomposition=(1, 1, comm.size, 1, 1))
with backend.open(context) as session:
    with session.run(calculation, scf=scf, material_id="Si",
                     restart=checkpoint, repartition=True) as state:
        assert state.scf_status["executed_iterations"] == 0
        state.continue_scf(steps=5)
```

New manifests record each file's three interior `grid_bounds` pairs. Admission
checks boxes, overlap/gaps, physical-input/runtime identity and file hashes.
Fortran maps accepted density, potentials, wavefunction/CG interiors, Fourier
density and every density/potential mixing history by global grid coordinates.
Periodic halo values are copied from saved interiors, without wavefunction
evolution. Mixing age, residual history, occupations and scalar state remain saved
values. Compare accepted interiors; old scratch halos need not have been current.

The existing zero-SCF bootstrap prepares the new topology's static grid, transform
and communication resources. Dead FFT partial-transform and reusable XC scratch
use that new layout rather than incompatible old buffers; they contain no SCF
history. No initial guess, density evaluation, eigenstate solve or SCF replay
reconstructs accepted values. Publication follows complete interior coverage and
a global success vote. A failed donor read preserves the live state.

Measured controls cover real/complex1↔2, x/y/z changes,16³2↔4 and
2×2×1→1×1×2 groups. Accepted cut interiors are exact; five subsequent iterations
match source-initialized continuous8 at rtol1e-9/atol1e-10. These are execution
checks, not physical convergence certification. The sample uses
`--spatial-parallel --save-scf NEW_CHECKPOINT` followed by
`--spatial-parallel --restart CHECKPOINT --repartition` with a different rank count.
Older saves without grid_bounds still support same-layout restoration.
HSE spatial GS uses the explicit entry below. ordinary mixed migration uses the entry below;
changed-runtime portability remains unsupported; common protected restore still keeps the same layout.

### Fast-PR HSE spatial rank-count restart

Use immutable `hse_spatial=True` on both hosts, including one rank, to select the
same full-support CPU Gamma collective exchange algorithm. It excludes
hse_orbital. HSE-spatial compute protocol1 and spatial-repartition protocol2 are
required; physical inputs, runtime, assets, material ID and SCF settings must match.

```python
backend = NativeGSBackend(build_record, work_root,
    decomposition=(1, 1, comm.size, 1, 1), hse_spatial=True)
with backend.open(context) as session:
    with session.run(calculation, scf=scf, material_id="Si",
                     restart=checkpoint, repartition=True) as state:
        assert state.scf_status["executed_iterations"] == 0
        state.continue_scf(steps=5)
```

Fortran maps dense ACE factor rows and exchange source/action/work/spatial-source/
previous rows by saved global physical grid coordinates, using the x-fast local
row layout. Replicated gauge, occupations, phase and mixing histories are retained.
Communicators are rebound; no accepted ACE evaluation, eigenstate/density solve
or SCF replay reconstructs the saved state. Static new-layout resources and dead
FFT/XC scratch follow the ordinary spatial contract above. Source data stays intact.

Fresh-process HSE1↔2, pre-PBE-stage1↔2,16³2↔4 and2×2×1→1×1×2 layouts passed.
Physical cut interiors match exactly, bootstrap SCF work is zero, and five new
iterations match the source-initialized continuous8 at the unchanged
rtol1e-9/atol1e-10. Late donor failures preserve originals before retry/continuation.
These are execution checks, not physical convergence certification.

The existing sample uses `--xc hse06 --spatial-parallel` on both sides; save with
`--save-scf NEW_CHECKPOINT`, then use `--restart CHECKPOINT --repartition` on the
new rank count. Existing default k/orbital and legacy CLI behavior is preserved.
The profile is periodic scalar, orthogonal-cell, fully occupied integer spin pairs,
fixed-ion, CPU, full-support Gamma with dense ACE and untiled gauges. Packed or
distributed ACE, tiled gauges, sparse localization and multi-k layouts; mixed orbital/spatial support requires the protocol2 entry below.
Changed-runtime migration and HSE spatial GS→RT remain outside this entry.

### Ordinary-XC mixed k/orbital/spatial SCF rank-count restart

A runtime declaring `sp_gs_scf_m_repartition_protocol()==1` supports ordinary
PZ/Libxc-PZ/PBE pure-to-mixed and mixed-to-mixed layout changes through the same
explicit `run(restart=..., repartition=True)` entry. The five decomposition
values are `(k, orbital, x, y, z)` and their product must equal the new host MPI
group size. Keep physical grid/k set, electrons/bands, SCF/mixing, material ID,
assets and runtime identical. Do not derive the physical k mesh from rank count.

```python
# Example source layout: (2, 2, 2, 1, 1), eight ranks.
# Example target layout: (1, 2, 1, 2, 1), four ranks.
backend = NativeGSBackend(build_record, work_root,
                         decomposition=(1, 2, 1, 2, 1))
with backend.open(context) as session:
    with session.run(calculation, scf=scf, material_id="Si",
                     restart=checkpoint, repartition=True) as state:
        assert state.scf_status["executed_iterations"] == 0
        state.continue_scf(steps=5)
```

Fortran treats each saved rank as a five-dimensional interior box: global x/y/z,
orbital and k indices from grid_bounds/orbital_bounds/k_bounds. All wavefunction/
CG slices are copied into the new ownership, including periodic halo copies.
Density, potential, Fourier density and mixing histories are already replicated
over old orbital/k groups: exactly one replica per spatial domain is used,
without adding it again. Native staging checks separate wave and scalar coverage
and publishes only after collective success. Static new-layout resources and
dead FFT/XC scratch follow the spatial contract; no density/eigenstate/SCF replay
reconstructs the saved accepted state. Default same-layout restore is unchanged.

Controls cover real/complex2→8,8→4,4→8 with two atoms, mixed4→2 and4→4 layout
changes. Physical cut wave/density values are exact; bootstrap SCF work is zero;
five new iterations match source-initialized continuous8 at the existing
rtol1e-9/atol1e-10. A missing late donor preserves original arrays and source
files. Independent mixed-layout owners and previous-runtime normal calculations
are also checked. These are execution controls, not physical convergence.

The existing sample accepts `--decomposition K O X Y Z` on both saving and
restoring hosts; its explicit-layout path keeps a fixed2×2×2 physical k mesh.
Use `--save-scf NEW_CHECKPOINT`, then `--restart CHECKPOINT --repartition` with
the new MPI count/layout. Fast HSE uses the Gamma protocol2 entry below; changed-runtime portability and
common protected restore with changed topology remain outside this entry.

### Fast-PR HSE Gamma orbital × spatial SCF restart

Use `hse_spatial=True` for both source and target hosts. Compute protocol
`sp_gs_hse_spatial_protocol()==2` permits simultaneous orbital and Cartesian
spatial decomposition, e.g. `(1,2,2,1,1)` (four ranks). The physical k mesh stays
Gamma `(1,1,1)`, without symmetry reduction, and every orbital group owns at
least one fully occupied spin-pair band. This extends the dense, full-support,
CPU, scalar, orthogonal-cell profile; it does not enable multi-k spatial HSE.

With mixed-repartition protocol2, use the existing explicit
`run(restart=checkpoint, repartition=True)` on a new rank group. Fortran maps
wavefunctions/CG history, dense ACE columns, cached exchange source/action and
spatial source/previous buffers by global orbital index and physical grid row.
Replicated gauge/occupation and scalar mixing state are retained without summing
replicas; communication handles belong to the new host. Restoration loads the
accepted numerical state without SCF/density/exchange reconstruction.

The existing example accepts `--xc hse06 --decomposition 1 O X Y Z` on both
sides; this path uses Gamma independently of MPI size. Save with `--save-scf`,
then resume with `--restart ... --repartition` into a new output directory.
Packed/distributed ACE, tiled gauges, sparse support, changed-runtime migration
and orbital/spatial GS→RT remain separate capabilities. Do not treat restart
regressions as physical convergence certification.

Measured controls cover2↔4,4→8,4→4 with changed spatial axes, uneven
orbital6→4 and pre-PBE cuts2↔4. Saved physical wave/density interiors are exact,
bootstrap SCF work is zero, and five new iterations agree with source-initialized
continuous8 at unchanged rtol1e-9/atol1e-10. Late donor read failures preserve
original state/files. Independent mixed owners and previous-runtime pure orbital/
spatial HSE and ordinary-XC controls verify the retained arithmetic.
