# Automatic dependency build validation — 2026-09-25

Platform: Apple Silicon, AppleClang 17 C, GNU Fortran 15, CMake 4.0.3,
Open MPI 5.0.9. Tests below concern build integration and bounded numerical
regressions, not new physical spectra or convergence studies.

- Ordinary CPU configure with only `CMAKE_BUILD_TYPE=Release`: HSE enabled,
  installed OpenBLAS/FFTW/Libxc selected, serial executable built. Fresh HSE
  GS integration input completed normally. No manually supplied compiler or
  dependency path and no extra C/Fortran flags.
- Missing-library path simulated with `CMAKE_DISABLE_FIND_PACKAGE_PkgConfig=TRUE`
  and `CMAKE_IGNORE_PREFIX_PATH=/opt/homebrew`, plus `USE_MPI=ON USE_LIBXC=ON`:
  downloaded and built FFTW 3.3.10, Libxc 5.2.3 including the legacy Fortran
  interface, and Netlib BLAS/LAPACK 3.12.1; linked one MPI executable.
- Earlier Nk=2^3 one-SCF-iteration GS → RT smoke cases 420/421: all six preparation/run/verification
  stages passed with the automatically built dependencies.
- Eight native numerical tests passed using the automatically built Libxc.
- Thirteen bounded input/restart checks passed, including default/explicit
  Taylor settings, inverse-length unit conversion, nondefault screening,
  invalid inputs, screening mismatch on restart, and non-HSE odd-step restart.
- `USE_HSE=OFF` serial build with forced Netlib fallback also built and passed
  the non-HSE restart checks in the preceding harness.

The fresh GS test exposed an Apple Accelerate `ZDOTC` complex-return ABI mismatch
with GNU Fortran. Automatic detection now selects OpenBLAS or builds Netlib for
that compiler/platform combination; rerunning the GS and RT tests passed.
The Libxc 5.x path also required updating all HSE external parameters together
and passing the correctly sized point count to the legacy spin interfaces.

Cross-platform CI and accelerator configurations have not been executed here.
Dependency source archives are hash pinned; generated libraries are confined to
the build tree. Existing vendor-specific LAPACK flags remain supported.

The earlier smoke cases did not establish SCF convergence. The stronger Nk=4^3
convergence tests subsequently exposed a separate Netlib eigenvector problem;
see [the updated portability and numerical test record](hse-platforms.md).

## Fugaku automatic selection preparation — 2026-09-27

Executed locally on Apple Silicon / GNU Fortran15, not on Fugaku:

- `python3 testsuites/unit_cmake/test_platform.py`: platform selection with fake compiler paths, overrides, compatibility alias, Release/Debug defaults, MPI/ScaLAPACK defaults, repeated inclusion, and dependency toolchain forwarding passed. No fake compiler was invoked as a real compiler.
- Fresh `cmake -S . -B <new-build> -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON`, then `cmake --build <new-build> -j4`: HSE enabled by default; installed OpenBLAS, Libxc and FFTW detected; full binary built. No manual compiler/library paths or extra compiler flags.
- Fresh `USE_HSE=OFF` Release configure and complete build passed.
- The newly built MPI/HSE executable passed `hse_lapack_eigenvectors` and `test_direct_wf.py` (MPI2/4, density/current/energy, finite support, ACE/U reuse, half dt, measured FFT and flag rejection).

The canonical `platforms/fugaku.cmake` is selected before `project()` on a Linux host with both Fujitsu MPI cross compiler wrappers available and no compiler/toolchain overrides. Legacy `fujitsu-a64fx-ea.cmake` remains an alias. Both default to vendor ScaLAPACK with MPI, with explicit OFF options respected. Target dependency configuration is compile/link only. Fujitsu compilation and vendor-library linkage still require testing on Fugaku; no remote test was executed in this session.
