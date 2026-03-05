set(CMAKE_Fortran_COMPILER      "mpifort")
set(CMAKE_C_COMPILER            "mpicc")
set(OPENMP_FLAGS                "-mp")

set(General_Fortran_FLAGS       "-Wall -fstrict-aliasing")
set(General_C_FLAGS             "-Wall -alias=ansi")

set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g -traceback -Mbounds ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g -traceback -Mbounds ${General_C_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O3 ${General_C_FLAGS}")

set(USE_MPI_DEFAULT             ON)

########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Compiling for x86_64")
set(CMAKE_SYSTEM_PROCESSOR "openmp")
