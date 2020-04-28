### Arm HPC compiler
set(ARCH_FLAGS                  "-mcpu=generic -march=armv8.1-a+sve")
set(OPENMP_FLAGS                "-fopenmp")
set(LAPACK_VENDOR_FLAGS         "-armpl -fno-simdmath") # NOTE: armpl option enables simdmath library but SVE version does not implemented.
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
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Arm HPC compiler for SVE processor (based Marvell ThunderX2)")
set(CMAKE_SYSTEM_PROCESSOR "armsve")
