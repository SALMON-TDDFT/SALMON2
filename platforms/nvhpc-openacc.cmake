set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_C_COMPILER            "mpicc")
set(OPENMP_FLAGS                "-Mnoopenmp")

set(General_Fortran_FLAGS       "-Wall -fstrict-aliasing")
set(General_C_FLAGS             "-Wall -alias=ansi")

set(OpenACC_FLAGS               "-acc=strict -gpu=cc70,cc80,managed,ptxinfo -Minfo=accel -DUSE_OPENACC")

set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g -traceback ${General_Fortran_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g -traceback ${General_C_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${General_Fortran_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O3 ${General_C_FLAGS} ${OpenACC_FLAGS}")

set(USE_MPI_DEFAULT             ON)
set(USE_OPENACC                 ON)

########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Compiling for x86_64 + Nvidia GPU")
set(CMAKE_SYSTEM_PROCESSOR "openacc")
