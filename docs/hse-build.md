# Building HSE

A C/Fortran compiler and CMake are required. For a normal CPU build:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
```

HSE is enabled by default. Add `-DUSE_MPI=ON` to configure an MPI build;
CMake uses the installed MPI implementation. No additional HSE dependency paths,
separate library build commands, or compiler workaround flags are needed.

CMake first checks installed dependencies, including their required interfaces.
Otherwise the build downloads verified source archives and builds static libraries
inside the build directory with the selected compilers:

| Dependency | Automatic source version | Purpose |
| --- | --- | --- |
| Libxc | 5.2.3 | HSE semilocal functional; optional conventional Libxc interface |
| FFTW | 3.3.10 | CPU exchange FFTs |
| Netlib BLAS/LAPACK | 3.12.1 | Fallback when a compatible installed implementation is absent |

The first fallback build requires network access. Installed compatible libraries
allow offline builds. Optional `LIBXC_INSTALLDIR`, `FFTW_INSTALLDIR`, and ordinary
CMake search paths remain available for managed installations. `USE_LIBXC=ON`
requires a compatible legacy Fortran module as well as the C library; CMake builds
both from the same source when the installed module is missing or incompatible.

On Apple systems with GNU Fortran, automatic selection uses OpenBLAS or builds
Netlib, avoiding the incompatible legacy Accelerate complex-return ABI. Explicit
vendor/toolchain settings remain the responsibility of that toolchain. GNU
Fortran 10+ compatibility flags for SALMON's legacy MPI interfaces are automatic.
C and Fortran OpenMP flags are detected separately.

`USE_HSE=OFF` omits HSE dependencies (Libxc is still needed if `USE_LIBXC=ON`).
OpenACC builds default to HSE disabled; enabling CPU HSE together with OpenACC
is unsupported. Cross-compiling toolchains are forwarded to HSE dependency
builds, but this change has only been executed and tested on Apple Silicon.

Source archives retain their upstream licenses: Libxc MPL-2.0, FFTW GPL-2.0-or-later,
and Netlib LAPACK modified BSD. See `LICENSE.THIRD-PARTY` and upstream `COPYING`
files in the downloaded sources.

The executed build and regression matrix is recorded in [validation](hse-build-validation.md).
