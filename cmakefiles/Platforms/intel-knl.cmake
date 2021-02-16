### Intel Compiler for Knights-Landing
set(ARCH_FLAGS                  "-xMIC-AVX512")
set(SIMD_SET                    "avx512")
set(OPENMP_FLAGS                "-qopenmp")
set(LAPACK_VENDOR_FLAGS         "-mkl=parallel")
set(ScaLAPACK_VENDOR_FLAGS      "-mkl=cluster")
set(Fortran_PP_FLAGS            "-fpp")

set(CMAKE_Fortran_COMPILER      "mpiifort")
set(CMAKE_C_COMPILER            "mpiicc")

set(General_FLAGS               "-qopt-ra-region-strategy=block -ansi-alias -fno-alias -qopt-report=5 -qopt-report-phase=openmp,loop,vec")
#set(General_Fortran_FLAGS       "-nogen-interface -std03 -warn all -diag-disable 6477,7025 ${General_FLAGS}")
set(General_Fortran_FLAGS       "-traceback -nogen-interface -std03 -warn all -diag-disable 6477,7025 ${General_FLAGS}")
set(General_C_FLAGS             "-Wall -diag-disable=10388 -restrict ${General_FLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g ${General_C_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O3 ${General_C_FLAGS}")

set(USE_MPI_DEFAULT                         ON)
set(USE_OPT_EXPLICIT_VECTORIZATION_DEFAULT  ON)


########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for Intel Knights-Landing")
set(CMAKE_SYSTEM_PROCESSOR "knl")
