### GNU compilers on macOS (native arm64/x86_64), Homebrew gcc + gfortran.
### Use:  mkdir build && cd build && python3 ../configure.py --arch=gnu-macos \
###         --disable-mpi --disable-scalapack --disable-libxc && make -j1
### First build must be -j1 (the legacy Fortran-module dep scan races under -j);
### incremental rebuilds can use -j<N>.
set(OPENMP_FLAGS                "-fopenmp")
set(Fortran_PP_FLAGS            "-cpp")

set(CMAKE_Fortran_COMPILER      "gfortran")
set(CMAKE_C_COMPILER            "gcc-15")

# -fallow-argument-mismatch: gfortran 10+ demotes rank/type mismatches to a
#   warning (the MPI/legacy interfaces in this tree rely on it).
# -Wno-error=implicit-function-declaration: gcc 14+ promotes implicit
#   declarations to an error; the POSIX C shim defines _XOPEN_SOURCE which hides
#   snprintf's prototype on this libc, so demote it back to a warning (the
#   function is still provided by the system C library at link time).
set(General_Fortran_FLAGS       "-ffree-line-length-none -fallow-argument-mismatch -w")
set(General_C_FLAGS             "-Wno-error=implicit-function-declaration -w")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O0 -g ${General_C_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O2 ${General_C_FLAGS}")

set(USE_MPI_DEFAULT             OFF)
