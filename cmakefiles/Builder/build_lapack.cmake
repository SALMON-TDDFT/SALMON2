include(ExternalProject)
include(${CMAKE_CURRENT_LIST_DIR}/hse_external_options.cmake)

if (LAPACK_VENDOR_FLAGS)
  message(STATUS "Set vendor-specific LAPACK libraries: ${LAPACK_VENDOR_FLAGS}")
  set(EXTERNAL_FLAGS ${LAPACK_VENDOR_FLAGS} ${EXTERNAL_FLAGS})
else ()
  # Apple's legacy Accelerate complex BLAS return ABI is incompatible with
  # GNU Fortran (e.g. ZDOTC). Prefer OpenBLAS, then build Netlib if absent.
  if(APPLE AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND NOT DEFINED BLA_VENDOR)
    set(BLA_VENDOR OpenBLAS)
    set(_salmon_saved_prefix_path "${CMAKE_PREFIX_PATH}")
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(PC_SALMON_BLAS QUIET openblas)
      if(PC_SALMON_BLAS_FOUND)
        list(APPEND CMAKE_PREFIX_PATH "${PC_SALMON_BLAS_PREFIX}")
      endif()
    endif()
    # Also cover package-manager keg-only libraries without pkg-config.
    find_library(SALMON_APPLE_OPENBLAS NAMES openblas
      PATH_SUFFIXES opt/openblas/lib openblas/lib)
    if(SALMON_APPLE_OPENBLAS)
      get_filename_component(_openblas_libdir "${SALMON_APPLE_OPENBLAS}" DIRECTORY)
      get_filename_component(_openblas_prefix "${_openblas_libdir}" DIRECTORY)
      list(APPEND CMAKE_PREFIX_PATH "${_openblas_prefix}")
    endif()
    set(_salmon_default_blas_vendor TRUE)
  endif()
  find_package(LAPACK QUIET)
  if(_salmon_default_blas_vendor)
    unset(BLA_VENDOR)
    set(CMAKE_PREFIX_PATH "${_salmon_saved_prefix_path}")
  endif()

  if (LAPACK_FOUND)
    message(STATUS "LAPACK library found.")
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
  else ()
    # NOTE: LAPACK 3.7.0 and later version can't build by GCC 4.8.5, which RHEL7 provided compiler.
    set(LAPACK_VERSION "3.12.1")
    message(STATUS "Build Netlib LAPACK library version ${LAPACK_VERSION}")

    # GNU 15/AArch64 testing found incorrect ZLARF1L eigenvectors when loop
    # vectorization was enabled, despite correct eigenvalues. Keep the
    # workaround local to the fallback, leaving SALMON/vendor BLAS optimized.
    set(_lapack_fortran_flags "${CMAKE_Fortran_FLAGS}")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
       AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 15
       AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 16
       AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|AARCH64)$")
      string(APPEND _lapack_fortran_flags " -fno-tree-loop-vectorize")
    endif()

    # old URL "http://www.netlib.org/lapack/lapack-${LAPACK_VERSION}.tgz"
    ExternalProject_Add(lapack-project
      URL              "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v${LAPACK_VERSION}.tar.gz"
      URL_HASH         SHA256=2ca6407a001a474d4d4d35f3a61550156050c48016d949f0da0529c0aa052422
      PREFIX           "${CMAKE_BINARY_DIR}/lapack"
      LIST_SEPARATOR   "|"
      BUILD_BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/lib/liblapack.a" "${CMAKE_CURRENT_BINARY_DIR}/lib/libblas.a"
      CMAKE_ARGS       ${SALMON_EXTERNAL_CMAKE_ARGS} -D BUILD_TESTING=off
                       -D CMAKE_INSTALL_LIBDIR=lib -D CMAKE_POLICY_VERSION_MINIMUM=3.5
                       -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
                       -D CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -D CMAKE_Fortran_FLAGS=${_lapack_fortran_flags}
                       -D CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}
                       -D CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}
      STEP_TARGETS     install
      EXCLUDE_FROM_ALL on
    )

    add_library(lapack STATIC IMPORTED)
    add_library(blas   STATIC IMPORTED)
    set_target_properties(lapack PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib/liblapack.a)
    set_target_properties(blas   PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib/libblas.a)
    add_dependencies(lapack lapack-project-install)
    add_dependencies(blas   lapack-project-install)
    set(EXTERNAL_LIBS lapack blas ${EXTERNAL_LIBS})
  endif ()
endif ()
