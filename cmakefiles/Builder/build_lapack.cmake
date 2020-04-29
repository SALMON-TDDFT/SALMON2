include(ExternalProject)

if (LAPACK_VENDOR_FLAGS)
  message(STATUS "Set vendor-specific LAPACK libraries: ${LAPACK_VENDOR_FLAGS}")
  set(EXTERNAL_FLAGS ${LAPACK_VENDOR_FLAGS} ${EXTERNAL_FLAGS})
else ()
  find_package(LAPACK QUIET)

  if (LAPACK_FOUND)
    message(STATUS "LAPACK library found.")
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
  else ()
    # NOTE: LAPACK 3.7.0 and later version can't build by GCC 4.8.5, which RHEL7 provided compiler.
    set(LAPACK_VERSION "3.6.1")
    message(STATUS "Build Netlib LAPACK library version ${LAPACK_VERSION}")

    ExternalProject_Add(lapack-project
      URL              "http://www.netlib.org/lapack/lapack-${LAPACK_VERSION}.tgz"
      PREFIX           "${CMAKE_BINARY_DIR}/lapack"
      CMAKE_ARGS       -D BUILD_SHARED_LIBS=off -D BUILD_TESTING=off
                       -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
                       -D CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -D CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
                       -D CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}
                       -D CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}
      STEP_TARGETS     install
      EXCLUDE_FROM_ALL on
    )

    add_library(lapack STATIC IMPORTED)
    add_library(blas   STATIC IMPORTED)
    set_target_properties(lapack PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib64/liblapack.a)
    set_target_properties(blas   PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib64/libblas.a)
    add_dependencies(lapack lapack-project-install)
    add_dependencies(blas   lapack-project-install)
    set(EXTERNAL_LIBS lapack blas ${EXTERNAL_LIBS})
  endif ()
endif ()
