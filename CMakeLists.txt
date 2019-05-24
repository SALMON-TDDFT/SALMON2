# recommended version : 3.0.2
cmake_minimum_required(VERSION 2.8)
enable_testing()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakefiles/Platforms)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakefiles/Modules)
include(${CMAKE_SOURCE_DIR}/cmakefiles/misc.cmake)
include(${CMAKE_SOURCE_DIR}/cmakefiles/setup.cmake)
include(${CMAKE_SOURCE_DIR}/cmakefiles/get_git_info.cmake)

### Project settings
project(SALMON Fortran C)


# Bulid options
## Library
option_set(USE_MPI                 "Use MPI parallelization" OFF)
option_set(USE_SCALAPACK           "Use ScaLAPACK Library"   OFF)

## Optimization
option_set(REDUCE_FOR_MANYCORE     "Enable optimization of reduction for many-core processor" OFF)
option_set(CURRENT_PREPROCESSING   "Enable preprocessing of the current computation in RT"    OFF)
option_set(LARGE_BLOCKING          "Enable large blocking for accelerator"                    OFF)

## Optimization for stencil compitations
option_set(OPT_STENCIL             "Enable optimization of the stencil"                                                ON)
option_set(DOMAIN_IS_POW2          "Enable whether the opimization assumes that a stencil domain size is power of two" OFF)
option_set(ARRAY_PADDING           "Enable array padding for the stencil"                                              ON)
option_set(EXPLICIT_VEC            "Enable explicit vectorization for the stencil"                                     OFF)
option_set(LOOP_BLOCKING           "Enable loop blocking at the stencil"                                               OFF)
option_set(SW_PREFETCH             "Enable software prefetch in the explicit vec of the stencil"                       OFF)

## test: joint hpsi subroutine
option_set(HPSI_TEST               "Use joint hpsi subroutine (test)"   OFF)
set(HPSI_TEST_SUFFIX "")
if (HPSI_TEST)
  set(HPSI_TEST_SUFFIX "_test")
endif()

## Additional library options
option_set(USE_LIBXC                 "Use LibXC library"                      OFF)
option_set(INSTALL_REQUIRED_PACKAGES "Request required packages installation" OFF)


# cross-compile or native-compile
if (${CMAKE_CROSSCOMPILING})
  # (For supercomputer systems) cross compile mode, use toolchain file.
  message(STATUS "Target platform ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
else ()
  # (For PC clusters and PC) compiler environment detects automatically.
  message(STATUS "Target platform is Native (${CMAKE_HOST_SYSTEM})")
  include(cmakefiles/find_build_environments.cmake)
endif ()

message(STATUS "user-defined FFLAGS = ${CMAKE_Fortran_FLAGS}")
message(STATUS "user-defined CFLAGS = ${CMAKE_C_FLAGS}")

set(TARGET_NAME         "salmon${TARGET_SUFFIX}")
set(CMAKE_Fortran_FLAGS "${ARCH} ${OPENMP_FLAGS} ${Fortran_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS} ${CMAKE_Fortran_FLAGS}")
set(CMAKE_C_FLAGS       "${ARCH} ${OPENMP_FLAGS} ${C_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS} ${CMAKE_C_FLAGS}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR})


add_definitions_if(USE_MPI               -DSALMON_USE_MPI)
add_definitions_if(CURRENT_PREPROCESSING -DARTED_CURRENT_PREPROCESSING)
add_definitions_if(REDUCE_FOR_MANYCORE   -DARTED_REDUCE_FOR_MANYCORE)
add_definitions_if(LARGE_BLOCKING        -DARTED_LBLK)

if (OPT_STENCIL)
  add_definitions(-DARTED_STENCIL_OPTIMIZED)

  add_definitions_if(DOMAIN_IS_POW2  -DARTED_DOMAIN_POWER_OF_TWO)
  add_definitions_if(EXPLICIT_VEC    -DARTED_EXPLICIT_VECTORIZATION)
  add_definitions_if(SW_PREFETCH     -DARTED_ENABLE_SOFTWARE_PREFETCH)
  add_definitions_if(ARRAY_PADDING   -DARTED_STENCIL_PADDING)
  add_definitions_if(LOOP_BLOCKING   -DARTED_STENCIL_ENABLE_LOOP_BLOCKING)
else ()
  add_definitions(-DARTED_STENCIL_ORIGIN)
endif ()
add_definitions("${ADDITIONAL_MACRO}")


if (LAPACK_INSTALLDIR)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${LAPACK_INSTALLDIR})
endif ()

if (LIBXC_INSTALLDIR)
  set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${LIBXC_INSTALLDIR})
endif ()

if (INSTALL_REQUIRED_PACKAGES)
  include(${CMAKE_SOURCE_DIR}/cmakefiles/install_required_packages.cmake)
else ()
  include(${CMAKE_SOURCE_DIR}/cmakefiles/find_required_packages.cmake)
endif ()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTERNAL_FLAGS}")
set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} ${EXTERNAL_FLAGS}")


include(${CMAKE_SOURCE_DIR}/cmakefiles/workaround.cmake)


# Generate configure header
configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.h.in
  ${CMAKE_BINARY_DIR}/version.h
)

configure_file(
  ${CMAKE_SOURCE_DIR}/src/versionf.h.in
  ${CMAKE_BINARY_DIR}/versionf.h
)

# source directory
set(SALMON_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
add_subdirectory(${SALMON_SRC_DIR})

# test directory
set(SALMON_TEST_DIR ${CMAKE_SOURCE_DIR}/testsuites)
add_subdirectory(${SALMON_TEST_DIR})

message(STATUS "Enabled macro")
get_directory_property(MACRO_LOG COMPILE_DEFINITIONS)
foreach(_macro ${MACRO_LOG})
  message(STATUS "  ${_macro}")
endforeach()
