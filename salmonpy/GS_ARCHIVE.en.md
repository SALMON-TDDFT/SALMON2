# Portable GS array archive (GS/2)

`salmonpy.gs_archive.save_ground_state` writes caller-supplied arrays; `GroundStateArchive.open` provides verified metadata and lazy reads. This is separate from native `GroundStateResult` (`salmonpy.gs/1`). Native GS→GS/2 export and SCF convergence certification are unimplemented. GS/2→TDDFT initialization is available under the explicit Gamma converter conditions below.

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

Bundled `include/salmonpy_gs_stream.f90` supports little-endian IEEE real64 hosts. `gs_read_chunk(path,values,status)` and `gs_write_chunk(path,values,status)` operate on flat real64 buffers. Status0 success,1 format/I/O/unsupported host,3 allocation failure,4 nonfinite data. Rejected reads leave caller values intact; writes refuse overwriting existing files. Big-endian hosts are rejected; Fujitsu compilation is not yet verified. The module does not implement JSON, hashes or full archive publication.

Python C-order(1,block,nz,ny,nx) chunk bytes match Fortran(nx,ny,nz,block), x fastest. Supply a flat buffer of nx*ny*nz*block values. SALMON native wfn.bin has stream branches but real/complex, MPI layout and spin conventions must be handled explicitly. info.bin/occupation.bin use compiler-dependent sequential records. Native GS/1→GS/2 export and multi-k/MPI conversion remain unsupported. Do not pass this archive directly as a native restart directory; use the explicit Gamma adapter below.

Format basis: [Zarr core specification](https://zarr-specs.readthedocs.io/en/latest/v3/core/index.html). Runtime dependencies remain NumPy only; an independent Zarr reader is used for synthetic compatibility validation.

Writer inputs are converted with NumPy: non-complex128 wavefunctions may require a full input-sized complex128 copy. Reading one band still reads a full orbital chunk; max_bytes bounds returned decoded bytes, not those temporary buffers. The Fortran helper uses int8 raw stream I/O so runtime REAL64 byte-swap settings cannot change chunk endian order; IEEE real64 and 8-bit file storage are checked.

## Gamma-only native TDDFT initialization

Use `GammaRestartConverter.from_build` and `archive.tddft_initial_state(converter=..., executable=..., pseudopotentials=..., max_steps=...)`. See the complete example in [GS_ARCHIVE.md](GS_ARCHIVE.md). Requires exact k=0/weight1, periodic unpolarized fixed ions, serial CPU runtime, grid>=8 per axis, even electrons and existing DFT XC scope(pz/libxc_pz/libxc_pbe). Libxc needs a matching runtime. Native density/potentials are reconstructed from orbitals/occupations; archived density/eigenvalues are not native restart records. Caller certifies physical GS provenance/convergence.

Build explicitly with `python tools/build_gs_restart_converter.py RUNTIME_BUILD_JSON NEW_DIRECTORY`; the installed script is under salmonpy/_runtime_source/tools. Requires a trusted audited serial GNU runtime record, retained source/CMake Unix Makefiles build files. Actual compiler and Fortran target flags are reused, application objects audited before link, and unchanged extracted native wavefunction reader/writer procedures test synthetic byte parity. No implicit build, download, runtime replacement or solver source import.

The converter is paired to one runtime SHA256. The supplied `Pseudopotential` objects must match archive species/hashes; local channels are explicit caller settings. `converter.write_bundle` only creates a new self-contained RT bundle and never starts a solver. `new_state` starts TDDFT. Temporary real/imaginary streams are removed after conversion. Source archive and existing destinations are preserved.

GNU GFORTRAN_CONVERT_UNIT global settings are recorded and must remain unchanged at every later spawn/replay; unit-specific settings are rejected. Compiler/runtime profile is caller-trusted provenance, not cryptographic publisher authentication. MPI runtimes, other target hashes and non-GNU profiles are rejected in this first version. Synthetic orthonormal orbital tests cover actual runtime initialization and bitwise independent/restored currents; they do not establish SCF convergence or real-material accuracy.

Potential hashes bind bytes/species, not the original Hamiltonian: local_channel, pseudopotential/XC consistency and physical GS provenance remain caller responsibilities. Conversion temporarily needs real+imaginary streams totaling16*grid_points*bands bytes in addition to the native complex output of the same size; these streams are removed afterward. Hashing native wavefunctions is streamed in1MiB buffers.
