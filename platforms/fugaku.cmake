### Fujitsu MPI cross compilers on the Fugaku login host
#set(OPENMP_FLAGS                "-Kopenmp")
set(OPENMP_FLAGS                "-Kopenmp -Nfjomplib")     # use Fujitsu OpenMP library
set(LAPACK_VENDOR_FLAGS         "-SSL2BLAMP")
set(ScaLAPACK_VENDOR_FLAGS      "-SCALAPACK -SSL2BLAMP")
set(Fortran_PP_FLAGS            "-Cpp")

set(CMAKE_Fortran_COMPILER      "mpifrtpx")
set(CMAKE_C_COMPILER            "mpifccpx")

# HSE caches use allocation/reallocation on intrinsic assignment.
set(General_Fortran_FLAGS       "-Kocl -Nlst=t -Koptmsg=2 -Ncheck_std=03s -Nalloc_assign")
set(General_C_FLAGS             "-Kocl -Nlst=t -Koptmsg=2 -Xg -std=gnu99")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g ${General_C_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-Kfast ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-Kfast ${General_C_FLAGS}")

set(USE_MPI_DEFAULT ON)
if(DEFINED USE_MPI AND NOT USE_MPI)
  set(USE_SCALAPACK_DEFAULT OFF)
else()
  set(USE_SCALAPACK_DEFAULT ON)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

########
# Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for Fujitsu A64FX")
set(CMAKE_SYSTEM_PROCESSOR    "a64fx")
set(CMAKE_Fortran_COMPILER_ID "Fujitsu" CACHE STRING "Fujitsu Fortran cross-compiler" FORCE)
set(CMAKE_C_COMPILER_ID       "Fujitsu" CACHE STRING "Fujitsu C cross-compiler" FORCE)
set(CMAKE_Fortran_MODDIR_FLAG "-M ")
