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
