cmake_minimum_required(VERSION 3.14)

if ("${CMAKE_CURRENT_SOURCE_DIR}" MATCHES "${CMAKE_CURRENT_BINARY_DIR}")
  message(FATAL_ERROR
  [[in-source build does not support
    Example:
      $ mkdir ${SALMON_BUILD_DIR}
      $ cd ${SALMON_BUILD_DIR}
      $ ${SALMON_SOURCE_DIR}/configure.py]]
  )
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/platforms)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakefiles/Modules)
include(${CMAKE_SOURCE_DIR}/cmakefiles/misc.cmake)


### Project settings
enable_testing()
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD          99)
project(SALMON
        VERSION   2.2.0
        LANGUAGES Fortran C
        )

message(STATUS "user-defined FFLAGS  = ${CMAKE_Fortran_FLAGS}")
message(STATUS "user-defined CFLAGS  = ${CMAKE_C_FLAGS}")
message(STATUS "user-defined LDFLAGS = ${CMAKE_EXE_LINKER_FLAGS}")

include(cmakefiles/check_build_environments.cmake)

if(USE_OPENACC)
  if(USE_CUDA)
    option(CMAKE_CUDA_ARCHITECTURES "80") # Dummy variable
    enable_language(CUDA)
  endif()
endif()


# Bulid options
## For third party libraries
option_set(USE_MPI       "Use MPI parallelization" OFF)
option_set(USE_SCALAPACK "Use ScaLAPACK Library"   OFF)
option_set(USE_EIGENEXA  "Use EigenExa Library"    OFF)
option_set(USE_LIBXC     "Use Libxc library"       OFF)
option_set(USE_FFTW      "Use FFTW library"        OFF)

## Optimization for stencil compitations
option_set(USE_OPT_ARRAY_PADDING           "Enable array padding for the stencil"          ON)
option_set(USE_OPT_EXPLICIT_VECTORIZATION  "Enable explicit vectorization for the stencil" OFF)


# variables
set(TARGET_NAME "salmon")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR})

include(${CMAKE_SOURCE_DIR}/cmakefiles/check_compiler_features.cmake)


# for third party libraries
include(${CMAKE_SOURCE_DIR}/cmakefiles/build_required_packages.cmake)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTERNAL_FLAGS}")
set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} ${EXTERNAL_FLAGS}")


# Generate configure header
configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.h.in
  ${CMAKE_BINARY_DIR}/version.h
)

configure_file(
  ${CMAKE_SOURCE_DIR}/src/versionf.h.in
  ${CMAKE_BINARY_DIR}/versionf.h
)

configure_file(
  ${CMAKE_SOURCE_DIR}/src/config.h.in
  ${CMAKE_BINARY_DIR}/config.h
)


add_subdirectory(${CMAKE_SOURCE_DIR}/src)        # sources
add_subdirectory(${CMAKE_SOURCE_DIR}/testsuites) # tests



# debug log
message(STATUS "Target executable file  = ${TARGET_NAME}")
message(STATUS "CMake generated FFLAGS  = ${CMAKE_Fortran_FLAGS}")
message(STATUS "CMake generated CFLAGS  = ${CMAKE_C_FLAGS}")
message(STATUS "CMake generated LDFLAGS = ${CMAKE_EXE_LINKER_FLAGS} ${EXTERNAL_LIBS}")
