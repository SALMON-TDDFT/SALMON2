### Fujitsu Compiler, Fugaku supercomputer (early access version)
set(ARCH                        "")
set(SIMD_SET                    "")
#set(OPENMP_FLAGS                "-Kopenmp")
set(OPENMP_FLAGS                "-Kopenmp -Nfjomplib")     # use Fujitsu OpenMP library
set(LAPACK_FLAGS                "-SSL2BLAMP")
set(ScaLAPACK_FLAGS             "-SCALAPACK -SSL2BLAMP")
set(ADDITIONAL_MACRO            "")
set(ADDITIONAL_OPTIMIZE_FLAGS   "")

set(Fortran_FLAGS_General       "-Kocl -Cpp -Nlst=t -Koptmsg=2 -Ncheck_std=03s")
set(C_FLAGS_General             "-Kocl      -Nlst=t -Koptmsg=2 -Xg -std=gnu99")

set(CMAKE_Fortran_COMPILER      "mpifrtpx")
set(CMAKE_C_COMPILER            "mpifccpx")

set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
set(CMAKE_C_FLAGS_DEBUG         "${CMAKE_Fortran_FLAGS_DEBUG}")

#set(CMAKE_Fortran_FLAGS_RELEASE "-Kfast -Kloop_nofission")  # stop loop fissing
set(CMAKE_Fortran_FLAGS_RELEASE "-Kfast")
set(CMAKE_C_FLAGS_RELEASE       "${CMAKE_Fortran_FLAGS_RELEASE}")

set(USE_MPI_DEFAULT ON)

########
# Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for Fujitsu A64FX")
set(CMAKE_SYSTEM_PROCESSOR    "a64fx")
set(CMAKE_Fortran_COMPILER_ID "Fujitsu" CACHE STRING "Fujitsu Fortran cross-compiler" FORCE)
set(CMAKE_C_COMPILER_ID       "Fujitsu" CACHE STRING "Fujitsu C cross-compiler" FORCE)
set(CMAKE_Fortran_MODDIR_FLAG "-M ")
