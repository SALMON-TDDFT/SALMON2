# Automatic HSE build dependencies

User requirement: ordinary CMake builds must prepare and link required libraries
without users specifying HSE library paths or workaround compiler flags.

1. Default HSE on for CPU, preserve explicit OFF and accelerator default OFF.
2. Reuse compatible installed FFTW/Libxc; automatically build pinned, hashed
   sources into the build directory if missing. Share Libxc with USE_LIBXC.
3. Keep existing BLAS/LAPACK vendor/system selection and make its fallback
   install layout deterministic. Preserve toolchain/cross-compile inputs.
4. Correct Libxc parameter updates for pre-7 libraries (set all parameters
   together), and include GNU legacy MPI compatibility automatically.
5. Validate system and bundled dependency paths with clean configure/build,
   HSE regression tests, nondefault omega and optional Libxc coexistence.
6. Update minimal user build instructions and dependency-license references.
