set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_C_COMPILER            "mpicc")
set(OPENMP_FLAGS                "-Mnoopenmp")

set(General_Fortran_FLAGS       "-Wall -fstrict-aliasing")
set(General_C_FLAGS             "-Wall -alias=ansi")

# OpenACC and CUDA
set(OpenACC_FLAGS               "-acc=strict -gpu=cc70,cc80,managed,ptxinfo -Minfo=accel -DUSE_OPENACC -DUSE_CUDA")
set(CMAKE_CUDA_FLAGS            "${CUDA_NVCC_FLAGS} -std=c++14 -rdc=false")
foreach(arch 70 80)
	set(CMAKE_CUDA_FLAGS          "${CMAKE_CUDA_FLAGS} -gencode arch=compute_${arch},code=sm_${arch}")
endforeach()

set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g -traceback ${General_Fortran_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g -traceback ${General_C_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${General_Fortran_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O3 ${General_C_FLAGS} ${OpenACC_FLAGS}")

set(USE_MPI_DEFAULT             ON)
set(USE_OPENACC                 ON)
set(USE_CUDA                    ON)

########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Compiling for x86_64 + Nvidia GPU")
set(CMAKE_SYSTEM_PROCESSOR "openacc")
