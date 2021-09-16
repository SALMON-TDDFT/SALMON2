set(CMAKE_Fortran_COMPILER      "mpifort")
set(CMAKE_C_COMPILER            "mpicc")
set(OPENMP_FLAGS                "-openmp")

set(General_Fortran_FLAGS       "-Wall -fstrict-aliasing")
set(OpenACC_FLAGS               "-DUSE_OPENACC -acc -Minfo -ta=tesla,managed -DUSE_CUDA")
set(General_C_FLAGS             "-Wall -restrict -ansi-alias -fno-alias")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g ${General_Fortran_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g ${General_C_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${General_Fortran_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O3 ${General_Fortran_FLAGS} ${OpenACC_FLAGS}")

set(USE_MPI_DEFAULT             ON)
set(USE_OPENACC ON)

########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Compiling for x86_64 + Nvidia GPU")
set(CMAKE_SYSTEM_PROCESSOR "openacc")

# CUDA
option(CUDA_USE_STATIC_CUDA_RUNTIME OFF)

set(CMAKE_CUDA_FLAGS  ${CUDA_NVCC_FLAGS} "-std=c++14 -arch=sm_80 -rdc=false")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++14 -O3")
set(CUDA_HOST_COMPILER nvcc)
