# HSE build portability and Nk=4^3 tests

## Status

The implementation uses the existing CPU toolchain selection. The Fugaku and
Linux configurations are source-audited, not certified by runs on those hosts:
this development session has no Fujitsu compiler or Linux runtime available.
Apple Silicon / GNU Fortran 15 is the executed platform. Passing its tests does
not replace tests with Fujitsu, Intel or a Linux GNU toolchain.

## Fugaku

Use the repository's existing toolchain rather than host compiler detection:

```sh
cmake -S . -B build-fugaku -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE="$PWD/platforms/fujitsu-a64fx-ea.cmake"
cmake --build build-fugaku -j 8
```

It selects mpifrtpx/mpifccpx, Fujitsu OpenMP, Linux/A64FX cross compilation and
SSL II (`-SSL2BLAMP`). The SSL II/OpenMP combination matches the
[RIKEN/RIST usage guide](https://www.r-ccs.riken.jp/fugaku/docs/workshop/2024/ja/introduction_to_fugaku_usage_seminar_ja_202410.pdf).
Those vendor BLAS flags bypass the Netlib fallback.
The dependency projects receive the selected compilers and the absolute toolchain
file. GNU-only compatibility/workaround flags are guarded by compiler identity.
The HSE library probes compile/link; they do not execute target programs on the
login host. Libxc native-host optimization is disabled, and FFTW is built serially
because SALMON calls its FFTs outside OpenMP regions.

Build-time source downloads require network access if compatible dependencies
are not already installed. A compute node is not required for these HSE library
probes, but MPI execution/CTest belongs inside the site's compute allocation.
Use the configured site's MPI launcher, four ranks and OMP_NUM_THREADS=12 for a
48-core node. That hybrid configuration still needs execution on Fugaku.

## Linux

For GNU C/Fortran plus installed MPI, ordinary CMake with `USE_MPI=ON` selects
compatible installed libraries or builds missing dependencies. For Intel oneAPI,
use `platforms/intel-oneapi.cmake` to retain its MPI/MKL/OpenMP settings. No
Apple-only OpenBLAS selection is applied on Linux. Linux x86_64 and AArch64
compiler/runtime tests remain required before claiming support verified there.

## Numerical dependency regression

During the Nk=4^3 test, the auto-built Netlib 3.12.1 path produced wrong ZHEEV
eigenvectors despite correct eigenvalues and INFO=0. A saved 16x16 SCF matrix
reproduced a normalized residual of 0.155; recompiling ZLARF1L without loop
vectorization reduced it to 7.1e-16. OpenBLAS gave 5.6e-16. An independent generated
Hermitian matrix reproduces the failure, without storing the SCF checkpoint.

The build now disables loop vectorization only for the Netlib fallback with
GNU Fortran 15.x on arm64/aarch64. This is a conservative architecture/compiler
guard based on the observed Mac failure; it is not a claim that the same failure
has been reproduced on Linux. The installed vendor-library and SALMON build flags
are unchanged. `hse_lapack_eigenvectors` checks the actual selected library's
residual and orthonormality before the GS fixture starts.

## Si integration tests

Cases 420/421 use a fresh HSE06 Si8 GS, 12^3 real-space grid, full shifted Nk=4^3,
16 occupied states, and a z-polarized impulse of 1e-4 au. GS convergence below
1e-8 is mandatory. Successful GS verification is the fixture required by RT
preparation, which also checks the producer's checkpoint and convergence log.

The 64-step Taylor4+ACE trajectory (dt=.16 au, 10.24 au total) is a CI regression.
It checks finite current/energy, bounded current, post-kick energy drift, and
the generated dielectric/conductivity relation. It does not resolve optical
peaks or establish k-point convergence.

```sh
OMP_NUM_THREADS=1 ctest --test-dir build \
  -R '(420_bulk_Si_hse_gs|421_bulk_Si_hse_rt)' --output-on-failure
```

Selecting only case 421 pulls in the LAPACK and GS fixtures automatically.
Python 3 is needed for verification, but not for building SALMON.

## Executed validation (2026-09-25/26)

Apple Silicon, GNU Fortran 15, AppleClang 17, Open MPI 5.0.9; Release builds.

| Dependency path | MPI ranks × OMP threads | GS | 64-step impulse |
| --- | --- | --- | --- |
| Automatically built Netlib 3.12.1 (workaround), Libxc 5.2.3, FFTW 3.3.10 | 4 × 1 | 110.47 s; residual 3.5022627e-9 | 170.04 s; verification passed |
| Installed OpenBLAS, Libxc 7, FFTW | 4 × 2 | 59.57 s; residual 3.4765297e-9 | 88.46 s; verification passed |

Both GS calculations converged at reported iteration 73. The installed-library
run selected only `421_bulk_Si_hse_rt`; CTest automatically executed all seven
LAPACK/GS/RT checks, with zero failures (150.37 s total). It used
`OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1`; SALMON reported two OpenMP threads.
The bundled run verified GS and RT in separate CTest invocations after correcting
a verifier that counted the trailing k-grid coefficient row as a k point.
These timings involve different libraries and thread counts, so they are not
an isolated OpenMP speedup measurement.

The GS/RT checks cover convergence, 64 k points with normalized weights, finite
outputs, the default HSE screening and Taylor propagator, electron count 32,
nonzero bounded impulse current, post-kick energy drift, and the dielectric/
conductivity identity. The original failing Netlib eigenvector path was reproduced
before applying the workaround. No new Fugaku or Linux execution is claimed.
