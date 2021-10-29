### Fujitsu Compiler, FX100 system
set(ARCH_FLAGS                  "-KHPC_ACE2")
set(OPENMP_FLAGS                "-Kopenmp")
set(LAPACK_VENDOR_FLAGS         "-SSL2BLAMP")
set(ScaLAPACK_VENDOR_FLAGS      "-SCALAPACK -SSL2BLAMP")
set(Fortran_PP_FLAGS            "-Cpp")

set(CMAKE_Fortran_COMPILER      "mpifrtpx")
set(CMAKE_C_COMPILER            "mpifccpx")

set(General_Fortran_FLAGS       "-Kocl,nooptmsg -v03s")
set(General_C_FLAGS             "-Kocl,nooptmsg -Xg -std=gnu99")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g ${General_C_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -Kfast,simd=1 ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O3 -Kfast,simd=1 ${General_C_FLAGS}")

set(USE_MPI_DEFAULT             ON)

########
# Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for Fujitsu SPARC64 XIfx")
set(CMAKE_SYSTEM_PROCESSOR    "s64fx")
set(CMAKE_Fortran_COMPILER_ID "Fujitsu" CACHE STRING "Fujitsu MPI Fortran cross-compiler" FORCE)
set(CMAKE_C_COMPILER_ID       "Fujitsu" CACHE STRING "Fujitsu MPI C cross-compiler" FORCE)
set(CMAKE_Fortran_MODDIR_FLAG "-M ")
