### Arm HPC compiler
set(TARGET_SUFFIX               ".cpu")

set(ARCH                        "-mcpu=generic -march=armv8.1-a+sve")
set(SIMD_SET                    "")
set(OPENMP_FLAGS                "-fopenmp")
set(LAPACK_FLAGS                "-armpl -fno-simdmath") # NOTE: armpl option enables simdmath library but SVE version does not implemented.
set(ScaLAPACK_FLAGS             "")
set(ADDITIONAL_MACRO            "-D__ARM_FLANG")
set(ADDITIONAL_OPTIMIZE_FLAGS   "-Wall -fstrict-aliasing")

set(Fortran_FLAGS_General       "-cpp")
set(C_FLAGS_General             "")

set(CMAKE_Fortran_COMPILER      "armflang")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffp-contract=fast -fma")
set(CMAKE_C_COMPILER            "armclang")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g")
set(CMAKE_C_FLAGS_RELEASE       "-O3 -ffp-contract=fast")

set(USE_MPI_DEFAULT             OFF)

########
# Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Arm HPC compiler for SVE processor (based Marvell ThunderX2)")
set(CMAKE_SYSTEM_PROCESSOR "armsve")
