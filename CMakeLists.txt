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

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakefiles/Platforms)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakefiles/Modules)
include(${CMAKE_SOURCE_DIR}/cmakefiles/misc.cmake)
include(${CMAKE_SOURCE_DIR}/cmakefiles/get_git_info.cmake)


### Project settings
enable_testing()
project(SALMON
        VERSION   2.0.0
        LANGUAGES Fortran C
        )

# cross-compile or native-compile
if (${CMAKE_CROSSCOMPILING})
  # (For supercomputer systems) cross compile mode, use toolchain file.
  message(STATUS "Target platform ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
else ()
  # (For PC clusters and PC) compiler environment detects automatically.
  message(STATUS "Target platform is Native (${CMAKE_HOST_SYSTEM})")
  include(cmakefiles/find_build_environments.cmake)
endif ()

message(STATUS "user-defined FFLAGS  = ${CMAKE_Fortran_FLAGS}")
message(STATUS "user-defined CFLAGS  = ${CMAKE_C_FLAGS}")
message(STATUS "user-defined LDFLAGS = ${CMAKE_EXE_LINKER_FLAGS}")


# Bulid options
## Library
option_set(USE_MPI                 "Use MPI parallelization" OFF)
option_set(USE_SCALAPACK           "Use ScaLAPACK Library"   OFF)

## Optimization
option_set(REDUCE_FOR_MANYCORE     "Enable optimization of reduction for many-core processor" OFF)
option_set(CURRENT_PREPROCESSING   "Enable preprocessing of the current computation in RT"    OFF)
option_set(LARGE_BLOCKING          "Enable large blocking for accelerator"                    OFF)

## Optimization for stencil compitations
option_set(DOMAIN_IS_POW2          "Enable whether the opimization assumes that a stencil domain size is power of two" OFF)
option_set(ARRAY_PADDING           "Enable array padding for the stencil"                                              ON)
option_set(EXPLICIT_VEC            "Enable explicit vectorization for the stencil"                                     OFF)
option_set(SW_PREFETCH             "Enable software prefetch in the explicit vec of the stencil"                       OFF)

## For third party libraries
option_set(USE_LIBXC               "Use Libxc library"                               OFF)
option_set(BUILD_REQUIRED_PACKAGES "Request build and static link required packages" OFF)

## test: joint hpsi subroutine
option_set(HPSI_TEST               "Use joint hpsi subroutine (test)"   OFF)
set(HPSI_TEST_SUFFIX "")
if (HPSI_TEST)
  set(HPSI_TEST_SUFFIX "_test")
endif()


# variables
set(TARGET_NAME         "salmon${TARGET_SUFFIX}") # FIXME: TARGET_SUFFIX will remove at version 2 ?
set(CMAKE_Fortran_FLAGS "${ARCH} ${OPENMP_FLAGS} ${Fortran_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS} ${CMAKE_Fortran_FLAGS}")
set(CMAKE_C_FLAGS       "${ARCH} ${OPENMP_FLAGS} ${C_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS} ${CMAKE_C_FLAGS}")
add_definitions("${ADDITIONAL_MACRO}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR})

add_compile_definitions_if(USE_MPI               SALMON_USE_MPI)

add_compile_definitions_if(EXPLICIT_VEC          SALMON_EXPLICIT_VECTORIZATION)
add_compile_definitions_if(SW_PREFETCH           SALMON_ENABLE_SOFTWARE_PREFETCH)
add_compile_definitions_if(DOMAIN_IS_POW2        SALMON_DOMAIN_POWER_OF_TWO)
add_compile_definitions_if(ARRAY_PADDING         SALMON_STENCIL_PADDING)

add_compile_definitions_if(CURRENT_PREPROCESSING ARTED_CURRENT_PREPROCESSING)
add_compile_definitions_if(REDUCE_FOR_MANYCORE   ARTED_REDUCE_FOR_MANYCORE)
add_compile_definitions_if(LARGE_BLOCKING        ARTED_LBLK)

include (${CMAKE_SOURCE_DIR}/cmakefiles/check_compiler_features.cmake)


# for third party libraries
if (LAPACK_INSTALLDIR)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${LAPACK_INSTALLDIR})
endif ()

if (LIBXC_INSTALLDIR)
  set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${LIBXC_INSTALLDIR})
endif ()

if (BUILD_REQUIRED_PACKAGES)
  include(${CMAKE_SOURCE_DIR}/cmakefiles/build_required_packages.cmake)
else ()
  include(${CMAKE_SOURCE_DIR}/cmakefiles/find_required_packages.cmake)
endif ()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTERNAL_FLAGS}")
set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} ${EXTERNAL_FLAGS}")

if (USE_MPI)
  find_package(MPI REQUIRED)
  if ("${MPI_Fortran_VERSION_MAJOR}" VERSION_LESS "3")
    message(FATAL_ERROR "SALMON requires MPI version 3 for nonblocking collectives.")
  endif ()
endif ()


# Generate configure header
configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.h.in
  ${CMAKE_BINARY_DIR}/version.h
)

configure_file(
  ${CMAKE_SOURCE_DIR}/src/versionf.h.in
  ${CMAKE_BINARY_DIR}/versionf.h
)


add_subdirectory(${CMAKE_SOURCE_DIR}/src)        # sources
add_subdirectory(${CMAKE_SOURCE_DIR}/testsuites) # tests


# debug log
message(STATUS "<<<BEGIN: BUILD ENVIRONMENT SUMMARY>>>")
message(STATUS "Target executable file  = ${TARGET_NAME}")
message(STATUS "CMake generated FFLAGS  = ${CMAKE_Fortran_FLAGS}")
message(STATUS "CMake generated CFLAGS  = ${CMAKE_C_FLAGS}")
message(STATUS "CMake generated LDFLAGS = ${CMAKE_EXE_LINKER_FLAGS} ${EXTERNAL_LIBS}")

message(STATUS "Enabled macro:")
get_directory_property(MACRO_LOG COMPILE_DEFINITIONS)
foreach(_macro ${MACRO_LOG})
  message(STATUS "  ${_macro}")
endforeach()
message(STATUS "<<<END: BUILD ENVIRONMENT SUMMARY>>>")
