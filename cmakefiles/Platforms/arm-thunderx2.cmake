### Arm HPC compiler
# NOTE: SALMON outputs NaN when compiling with Arm HPC compiler version 19.1, don't use it.
set(TARGET_SUFFIX               ".cpu")

set(ARCH                        "-mcpu=thunderx2t99")
set(SIMD_SET                    "")
set(OPENMP_FLAGS                "-fopenmp")
set(LAPACK_FLAGS                "-armpl")
set(ScaLAPACK_FLAGS             "")
set(ADDITIONAL_MACRO            "-D__ARM_FLANG")
set(ADDITIONAL_OPTIMIZE_FLAGS   "-Wall -fstrict-aliasing")

set(Fortran_FLAGS_General       "-cpp")
set(C_FLAGS_General             "")

set(CMAKE_Fortran_COMPILER      "armflang")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -ffp-contract=fast")
set(CMAKE_C_COMPILER            "armclang")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g")
set(CMAKE_C_FLAGS_RELEASE       "-Ofast -ffp-contract=fast")

set(USE_MPI             OFF)
set(REDUCE_FOR_MANYCORE ON)

########
# Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Arm HPC compiler for Marvell ThunderX2")
set(CMAKE_SYSTEM_PROCESSOR "thunderx2")
