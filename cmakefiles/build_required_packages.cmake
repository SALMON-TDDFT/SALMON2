include(ExternalProject)


### Numerical library
if (USE_SLACAPACK)
  message(FATAL_ERROR "We don't support builiding ScaLAPACK. Please contact your system administrator.")
else ()
  set(FILE_MATHLIB lapack)
  find_package(LAPACK QUIET)

  if (LAPACK_FOUND)
    message(STATUS "LAPACK library found.")
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${LAPACK_LIBRARIES})
    set(EXTERNAL_FLAGS ${EXTERNAL_FLAGS} ${LAPACK_FLAGS})
  else ()
    set(LAPACK_VERSION "3.8.0")
    message(STATUS "Build Netlib LAPACK library version ${LAPACK_VERSION}")

    ExternalProject_Add(lapack-project
      URL              "http://www.netlib.org/lapack/lapack-${LAPACK_VERSION}.tar.gz"
      PREFIX           "${CMAKE_BINARY_DIR}/lapack"
      CMAKE_ARGS       -D BUILD_SHARED_LIBS=off
                       -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
                       -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER} -D CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -D CMAKE_C_FLAGS=${CMAKE_C_FLAGS} -D CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
                       -D CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG} -D CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}
                       -D CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE} -D CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}
      STEP_TARGETS     install
      EXCLUDE_FROM_ALL on
    )

    add_library(lapack STATIC IMPORTED)
    add_library(blas   STATIC IMPORTED)
    set_target_properties(lapack PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib64/liblapack.a)
    set_target_properties(blas   PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib64/libblas.a)
    add_dependencies(lapack lapack-project-install)
    add_dependencies(blas   lapack-project-install)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} lapack blas)
  endif ()
endif ()


### Libxc
if (USE_LIBXC)
  find_package(Libxc QUIET)

  if (Libxc_FOUND)
    message(STATUS "Libxc found.")
    include_directories($<TARGET_PROPERTY:Libxc::xc,INTERFACE_INCLUDE_DIRECTORIES>)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} Libxc::xc)
  else ()
    set(LIBXC_VERSION "4.3.4")
    message(STATUS "Build Libxc version ${LIBXC_VERSION}")

    ExternalProject_Add(libxc-project
      URL              "http://www.tddft.org/programs/octopus/down.php?file=libxc/${LIBXC_VERSION}/libxc-${LIBXC_VERSION}.tar.gz"
      PREFIX           "${CMAKE_BINARY_DIR}/libxc"
      CMAKE_ARGS       -D BUILD_SHARED_LIBS=off -D ENABLE_FORTRAN=on -D ENABLE_XHOST=off
                       -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
                       -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER} -D CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -D CMAKE_C_FLAGS=${CMAKE_C_FLAGS} -D CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
                       -D CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG} -D CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}
                       -D CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE} -D CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}
      STEP_TARGETS     install
      EXCLUDE_FROM_ALL on
    )

    add_definitions(-DSALMON_USE_LIBXC)
    add_library(xcf90 STATIC IMPORTED)
    add_library(xc    STATIC IMPORTED)
    set_target_properties(xcf90 PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib64/libxcf90.a)
    set_target_properties(xc    PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib64/libxc.a)
    add_dependencies(xcf90 libxc-project-install)
    add_dependencies(xc    libxc-project-install)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} xcf90 xc)
  endif ()
endif ()


# add search path
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
link_directories("${CMAKE_CURRENT_BINARY_DIR}/lib")
