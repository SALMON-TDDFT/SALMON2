### Arm HPC compiler
# NOTE: SALMON outputs NaN when compiling with Arm HPC compiler version 19.1, don't use it.
set(ARCH_FLAGS                  "-mcpu=thunderx2t99")
set(OPENMP_FLAGS                "-fopenmp")
set(LAPACK_VENDOR_FLAGS         "-armpl")
set(Fortran_PP_FLAGS            "-cpp")

set(CMAKE_Fortran_COMPILER      "armflang")
set(CMAKE_C_COMPILER            "armclang")

set(General_Fortran_FLAGS       "-Wall -fstrict-aliasing")
set(General_C_FLAGS             "${General_Fortran_FLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g ${General_C_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffp-contract=fast -fma ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O3 -ffp-contract=fast ${General_C_FLAGS}")

########
# Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Arm HPC compiler for Marvell ThunderX2")
set(CMAKE_SYSTEM_PROCESSOR "thunderx2")
