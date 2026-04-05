cmake_minimum_required(VERSION 3.14)

if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR
  [[in-source build does not support
    Example:
      $ mkdir ${SALMON_BUILD_DIR}
      $ cd ${SALMON_BUILD_DIR}
      $ ${SALMON_SOURCE_DIR}/configure.py]]
  )
endif ()

include(${CMAKE_SOURCE_DIR}/cmakefiles/misc.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/platforms)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakefiles/Modules)

### Project settings
enable_testing()
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD          99)
project(SALMON
        VERSION   2.2.2
        LANGUAGES Fortran C
        )

if(DEFINED CMAKE_Fortran_FLAGS_DEBUG)
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" CACHE STRING "Fortran debug flags" FORCE)
endif()
if(DEFINED CMAKE_Fortran_FLAGS_RELEASE)
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" CACHE STRING "Fortran release flags" FORCE)
endif()
if(DEFINED CMAKE_C_FLAGS_DEBUG)
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}" CACHE STRING "C debug flags" FORCE)
endif()
if(DEFINED CMAKE_C_FLAGS_RELEASE)
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "C release flags" FORCE)
endif()

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
option(USE_MPI       "Use MPI parallelization" OFF)
option(USE_SCALAPACK "Use ScaLAPACK Library"   OFF)
option(USE_EIGENEXA  "Use EigenExa Library"    OFF)
option(USE_LIBXC     "Use Libxc library"       OFF)
option(USE_FFTW      "Use FFTW library"        OFF)

## Optimization for stencil compitations
option(USE_OPT_ARRAY_PADDING           "Enable array padding for the stencil"          ON)
option(USE_OPT_EXPLICIT_VECTORIZATION  "Enable explicit vectorization for the stencil" OFF)


# variables
set(TARGET_NAME "salmon")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR})

# Enable Fortran preprocessor for free form sources
#set(CMAKE_Fortran_PREPROCESS_SOURCE ON)
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffree-form -x f95-cpp-input")

# Respect platform-specific preprocessing flags when provided by the toolchain.
if(DEFINED Fortran_PP_FLAGS)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Fortran_PP_FLAGS}")
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
endif()

include(${CMAKE_SOURCE_DIR}/cmakefiles/check_compiler_features.cmake)


# for third party libraries
include(${CMAKE_SOURCE_DIR}/cmakefiles/build_required_packages.cmake)

# Find and add OpenMP
find_package(OpenMP REQUIRED Fortran)
find_package(OpenMP QUIET)
if(OpenMP_Fortran_FOUND)
  list(APPEND EXTERNAL_LIBS OpenMP::OpenMP_Fortran)
endif()
if(OpenMP_C_FOUND)
  list(APPEND EXTERNAL_LIBS OpenMP::OpenMP_C)
endif()

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


# Add subdirectory for sources
add_subdirectory(${CMAKE_SOURCE_DIR}/src)        # sources
add_subdirectory(${CMAKE_SOURCE_DIR}/testsuites) # tests




# debug log
message(STATUS "Target executable file  = ${TARGET_NAME}")
message(STATUS "CMake generated FFLAGS  = ${CMAKE_Fortran_FLAGS}")
message(STATUS "CMake generated CFLAGS  = ${CMAKE_C_FLAGS}")
message(STATUS "CMake generated LDFLAGS = ${CMAKE_EXE_LINKER_FLAGS} ${EXTERNAL_LIBS}")
