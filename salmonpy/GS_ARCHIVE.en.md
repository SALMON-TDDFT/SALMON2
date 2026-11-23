# Portable GS array archive (GS/2)

`salmonpy.gs_archive.save_ground_state` writes caller-supplied arrays; `GroundStateArchive.open` provides verified metadata and lazy reads. This is separate from solver `GroundStateResult` (`salmonpy.gs/1`). Explicit export of DFT-generated Gamma GS/1 is supported below; SCF is not rerun or recertified. GS/2→TDDFT initialization is available under the explicit Gamma converter conditions below.

```python
from salmonpy.gs_archive import GroundStateArchive
archive = GroundStateArchive.open('existing-gs-array-archive')
u = archive.wavefunctions(0, 0, 4)  # k=0, bands[0,4), shape(band,z,y,x)
occupations = archive.array('occupations')
```

See the synthetic write example in [GS_ARCHIVE.md](GS_ARCHIVE.md). It demonstrates interchange, not a converged material. Do not mutate inputs during writing or edit an archive while reading. Only new destinations are accepted; failed writes remove their in-progress directory. The manifest is written last so incomplete archives cannot open.

## Physical conventions

Initial scope: periodic, unpolarized real-space arrays. `wavefunctions[k,band,z,y,x]` is complex128 cell-periodic u, with `psi=exp(+i*k_cartesian dot r)*u`. Rows of `cell_bohr` are lattice vectors a1,a2,a3; reciprocal rows are `B=2*pi*inverse(A).T`. kpoints are coefficients of those reciprocal vectors, not Cartesian components. The real-space grid starts at fractional origin zero and excludes duplicated periodic endpoints:

`r(ix,iy,iz)=(ix/nx)*a1+(iy/ny)*a2+(iz/nz)*a3`, indices from0 through dimension-1.

Array/grid order is(z,y,x), so x varies fastest in C bytes. Atomic positions are fractional lattice coefficients. Units: cell bohr, u bohr^(-3/2), density bohr^(-3), eigenvalues/energy hartree. Nonnegative k weights sum to one; occupations lie in[0,2]; their weighted sum equals electrons. Density is nonnegative, cell volume positive. Finite values are checked. Wavefunction normalization/orthogonality, density consistency and SCF convergence are not certified.

## On-disk contract

Root `manifest.json` contains schema `salmonpy.gs/2`, metadata, exact layouts and SHA256 of every data/Zarr-metadata file. Root `zarr.json` is a Zarr-v3 group. Arrays: wavefunction_real, wavefunction_imag, density, eigenvalues, occupations, kpoints, weights, cell_bohr, fractional_positions, atomic_numbers. Atomic numbers are int64; others float64, all uncompressed bytes codec, little endian, C order. Separate real/imaginary arrays avoid complex-dtype extensions.

Wave chunks have shape(1,orbital_block,nz,ny,nx); eigenvalues/occupations use(1,orbital_block). Other arrays occupy one chunk. Edge chunks retain full size and zero-fill out-of-array values. Default orbital_block16 is an I/O granularity, not a memory limit, and is capped to the band count. No spatial partitioning in this version.

Pseudopotentials contain references only: atomic_number, filename, sha256, source. No potential bytes or automatic download. Multiple different potentials for the same atomic number are unsupported.

`open` validates every file's layout/hash/finiteness without allocating the full wavefunction. Metadata and density arrays plus chunk buffers are used. `wavefunctions(k,start,stop)` reads only relevant blocks and rechecks hashes; returned metadata/arrays are immutable. Unexpected files/directories, symlinks, unsupported layouts and stale chunks are rejected.

No arbitrary machine-wide capacity ceiling. Optional `open(max_bytes=...)` bounds total archive-file bytes; array/wavefunction max_bytes bounds returned decoded bytes. Neither bounds process RSS: chunk buffers, validation and immutable copies add storage. Hashes detect corruption, not publisher authenticity.

## Fortran stream interface and limits

Bundled `include/salmonpy_gs_stream.f90` supports little-endian IEEE real64 hosts. `gs_read_chunk(path,values,status)` and `gs_write_chunk(path,values,status)` operate on flat real64 buffers. Status0 success,1 format/I/O/unsupported calling program,3 allocation failure,4 nonfinite data. Rejected reads leave caller values intact; writes refuse overwriting existing files. Big-endian hosts are rejected; Fujitsu compilation is not yet verified. The module does not implement JSON, hashes or full archive publication.

Python C-order(1,block,nz,ny,nx) chunk bytes match Fortran(nx,ny,nz,block), x fastest. Supply a flat buffer of nx*ny*nz*block values. SALMON solver wfn.bin has stream branches but real/complex, MPI layout and spin conventions must be handled explicitly. info.bin/occupation.bin use compiler-dependent sequential records. solver GS/1→GS/2 export supports the strict Gamma conditions below; multi-k/MPI conversion remains unsupported. Do not pass this archive directly as a solver restart directory; use the explicit Gamma adapter below.

Format basis: [Zarr core specification](https://zarr-specs.readthedocs.io/en/latest/v3/core/index.html). Runtime dependencies remain NumPy only; an independent Zarr reader is used for synthetic compatibility validation.

Writer inputs are converted with NumPy: non-complex128 wavefunctions may require a full input-sized complex128 copy. Reading one band still reads a full orbital chunk; max_bytes bounds returned decoded bytes, not those temporary buffers. The Fortran helper uses int8 raw stream I/O so runtime REAL64 byte-swap settings cannot change chunk endian order; IEEE real64 and 8-bit file storage are checked.

## Gamma-only TDDFT called directly from Python initialization

Use `GammaRestartConverter.from_build` and `archive.tddft_initial_state(converter=..., executable=..., pseudopotentials=..., max_steps=...)`. See the complete example in [GS_ARCHIVE.md](GS_ARCHIVE.md). Requires exact k=0/weight1, periodic unpolarized fixed ions, serial CPU runtime, grid>=8 per axis, even electrons and existing DFT XC scope(pz/libxc_pz/libxc_pbe). Libxc needs a matching runtime. solver density/potentials are reconstructed from orbitals/occupations; archived density/eigenvalues are not solver restart records. Caller certifies physical GS provenance/convergence.

Build explicitly with `python tools/build_gs_restart_converter.py RUNTIME_BUILD_JSON NEW_DIRECTORY`; the installed script is under salmonpy/_runtime_source/tools. Requires a trusted audited serial GNU runtime record, retained source/CMake Unix Makefiles build files. Actual compiler and Fortran target flags are reused, application objects audited before link, and unchanged extracted solver wavefunction reader/writer procedures test synthetic byte parity. No implicit build, download, runtime replacement or solver source import.

The converter is paired to one runtime SHA256. The supplied `Pseudopotential` objects must match archive species/hashes; local channels are explicit caller settings. `converter.write_bundle` only creates a new self-contained RT bundle and never starts a solver. `new_state` starts TDDFT. Temporary real/imaginary streams are removed after conversion. Source archive and existing destinations are preserved.

GNU GFORTRAN_CONVERT_UNIT global settings are recorded and must remain unchanged at every later spawn/replay; unit-specific settings are rejected. Compiler/runtime configuration is caller-trusted provenance, not cryptographic publisher authentication. MPI runtimes, other target hashes and non-GNU configurations are rejected in this first version. Synthetic orthonormal orbital tests cover actual runtime initialization and bitwise independent/restored currents; they do not establish SCF convergence or real-material accuracy.

Potential hashes bind bytes/species, not the original Hamiltonian: local_channel, pseudopotential/XC consistency and physical GS provenance remain caller responsibilities. Conversion temporarily needs real+imaginary streams totaling16*grid_points*bands bytes in addition to the solver complex output of the same size; these streams are removed afterward. Hashing solver wavefunctions is streamed in1MiB buffers.

## Export existing Gamma GS/1

Use `gs = GroundStateResult.open('existing-native-gs-result')`, then `archive = gs.export_archive('new-portable-gs', converter=converter, orbital_block=16)`. Explicitly rebuild a new converter with solver export support; old import-only build records remain usable for import but cannot export. Existing outputs/binaries are not overwritten.

Only SALMONpy DFT-generated Γ/serial/unpolarized GS/1 with one real or complex wfn.bin is supported. Verify full source hashes, paired runtime, header k/band/iteration/rank counts, occupations/electrons, physical input groups and metadata. Hand-edited/foreign input layouts are rejected. The known wavefunction-only warm-start input variant is supported; mixing history is not migrated. Eigenvalues come from native/gs_eigen.data in atomic units; missing/eV/nonfinite/inconsistent tables are rejected. Table occupations must equal solver occupation records.

Density is reconstructed in band order as `sum_b occ[b]*abs(u_b(r))**2`, rather than copied from the mixed solver density. Reported energies/eigenvalues are retained, not recomputed. This is not a complete SCF-state checkpoint. Optional typed metadata.provenance records source/runtime/exporter hashes, derived-density origin, solver orbital kind, local channels and I/O environment. Existing GS/2 without provenance still loads; reimport enforces recorded local channels. Provenance/hashes do not authenticate a publisher or certify physical correctness.

Source GS is unchanged; three verified records are copied privately. solver export holds one orbital; Python maps split streams and a complex scratch file rather than heap-allocating all wavefunctions. OS-managed mapped pages can still consume RSS; temporary disk is required. Mappings/temp files are released after success/failure. Only a newly claimed incomplete output is removed on failure. Synthetic solver real/complex and end-to-end TDDFT tests validate value/current/replay parity, not converged-material SCF accuracy.

## Opt-in serial multi-k conversion

`SerialRestartConverter.from_build('NEW_BUILD/build.json')` enables multi-k support only with a newly built/probed multi-k-capable GNU converter. Existing `GammaRestartConverter` remains exact-Gamma and old import-only records keep their previous behavior. Pass the new converter to the same `GroundStateResult.export_archive` and `GroundStateArchive.tddft_initial_state` methods.

Current multi-k scope: axis-aligned positive-diagonal cells, periodic unpolarized fixed ions, one serial CPU rank. Export supports the strict DFT-generated regular mesh, requires native/gs_k.data in its orthogonal reduced-coordinate format, verifies reciprocal vectors, point order/coordinates and uniform weights against input kmesh (k coordinate atol1e-14; weight rtol1e-14). Nonorthogonal/Cartesian k exports and custom solver file_kw source exports are rejected. Import supports caller-supplied reduced k coordinates and normalized nonnegative weights, including nonuniform distributions.

solver files keep k-major/band-minor order and occupations[band,k] as Fortran records. Portable arrays keep(k,band,z,y,x); density becomes `sum_k,b weights[k]*occupations[k,b]*abs(u[k,b])**2` in fixed k/band order. solver conversion checks weighted electron count and16*gridpoints*bands*kpoints sizes, with64-bit indices/capacity arithmetic. One orbital's solver scratch is retained; occupation/weight storage scales with k×bands.

Import creates an explicit reduced-coordinate `kpoints.dat` referenced by file_kw, whose hash is part of the conversion validation record and checked before loading. GNU record environment binding and source/archive protection remain in force. Weights are normalized by the solver file_kw loader. Exported GS/2 provenance adds optional kpoint_source='native/gs_k.data'.

Wavefunctions and occupations retain their solver bits. solver k-table text has limited decimal precision and file_kw renormalizes weights, so general k/weight/current bit equality is not promised. Actual-runtime synthetic tests use rtol1e-12,atol1e-14 for original vs reimport currents.2/8 grids and nonuniform2-point import matched bitwise;3/6 grids gave differences around1e-18. Restoring the same imported trajectory matched bitwise. These are synthetic integration tests, not SCF/real-material accuracy certification. MPI/spin/nonorthogonal/other compilers remain future work.

```python
from salmonpy.gs_native import SerialRestartConverter
converter = SerialRestartConverter.from_build('NEW_MULTIK_CONVERTER/build.json')
archive = gs.export_archive('new-portable-multik-gs', converter=converter)
# Same explicit paired engine and matching potentials as the Gamma example:
initial = archive.tddft_initial_state(converter=converter, executable=engine,
                                    pseudopotentials=potentials, max_steps=1000)
```

Multi-k states support the existing `state.save_checkpoint('new-history')` and `state.restore_checkpoint('new-history')` on a freshly created identical model with fixed dt. The GS/2 archive and paired converter/runtime/potentials must be supplied separately; checkpoints store history, without duplicating GS arrays. Synthetic 2/3/6/8-point and nonuniform two-point tests saved after six steps, restored in a fresh Python process, and matched the following 20 steps bitwise. Two states with different interleaved input histories, including continuation after closing a sibling, matched independent runs bitwise. Each state owns a separate worker process; these results do not establish multiple TDDFT called directly from Python contexts inside one process.

## Retained calculation objects

`NativeGSState.export_archive` adds a collective regular-mesh DFT-owned export entry (root writes/returns the archive).
See [NATIVE_GS.en.md](NATIVE_GS.en.md) for scope, provenance and measured limits.
It is not a complete SCF checkpoint or solver GS-handle restoration.

## solver retained-GS initialization

`NativeGSSession.load_archive` converts GS/2 once on subgroup root and imports a
new calculation object with no SCF. Explicit converter/runtime/potential selection is
required. See [NATIVE_GS.en.md](NATIVE_GS.en.md) and examples/portable_gs.py.
This restores array initialization, not archived energy/eigenvalues, mixing
history or an in-progress solver SCF checkpoint.

GS/2-imported retained owners can re-export using
`salmonpy.retained-gs-reexport/1` provenance. solver orbitals/occupations and
reconstructed density are separated from carried source energy/eigenvalues/model/
sampling information. Current convergence remains unknown. See the retained-GS
guide; existing provenance stays readable, and new provenance needs a matching reader.

TB-mBJ archives retain root metadata `xc_parameters: {"cval": null}` for automatic mode or a fixed coefficient. solver conversion rejects missing TB-mBJ parameters. Custom XC parameter dictionaries may contain finite JSON values; the built-in solver converter rejects unsupported parameters instead of silently dropping them.
