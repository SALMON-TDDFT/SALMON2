include(ExternalProject)


if (${CMAKE_CROSSCOMPILING})
  message(FATAL_ERROR "When cross compile mode, we don't support the automatic installation of required packages. Please contact your system administrator.")
endif ()


### Numerical library
if (USE_SLACAPACK)
  message(FATAL_ERROR "We don't support an automatic installation of ScaLAPACK. Please contact your system administrator.")
else ()
  set(FILE_MATHLIB lapack)
  find_package(LAPACK QUIET)

  if (LAPACK_FOUND)
    message(STATUS "LAPACK library found.")
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${LAPACK_LIBRARIES})
    set(EXTERNAL_FLAGS ${EXTERNAL_FLAGS} ${LAPACK_FLAGS})
  else ()
    set(LAPACK_VERSION "3.8.0")
    message(STATUS "Install Netlib LAPACK library version ${LAPACK_VERSION}")

    ExternalProject_Add(lapack-project
      URL              "http://www.netlib.org/lapack/lapack-${LAPACK_VERSION}.tar.gz"
      PREFIX           "${CMAKE_BINARY_DIR}/lapack"
      CMAKE_ARGS       -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      STEP_TARGETS     install
      EXCLUDE_FROM_ALL on
    )

    add_library(lapack STATIC IMPORTED)
    add_library(blas   STATIC IMPORTED)
    set_target_properties(lapack PROPERTIES IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib64/liblapack.a)
    set_target_properties(blas   PROPERTIES IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib64/libblas.a)
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
    message(STATUS "Enable installation Libxc version ${LIBXC_VERSION}")

    ExternalProject_Add(libxc-project
      URL              "http://www.tddft.org/programs/octopus/down.php?file=libxc/${LIBXC_VERSION}/libxc-${LIBXC_VERSION}.tar.gz"
      PREFIX           "${CMAKE_BINARY_DIR}/libxc"
      CMAKE_ARGS       -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -D ENABLE_FORTRAN=on
      STEP_TARGETS     install
      EXCLUDE_FROM_ALL on
    )

    add_definitions(-DSALMON_USE_LIBXC)
    add_library(xcf90 STATIC IMPORTED)
    add_library(xc    STATIC IMPORTED)
    set_target_properties(xcf90 PROPERTIES IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib64/libxcf90.a)
    set_target_properties(xc    PROPERTIES IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib64/libxc.a)
    add_dependencies(xcf90 libxc-project-install)
    add_dependencies(xc    libxc-project-install)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} xcf90 xc)
  endif ()
endif ()


# add search path
include_directories("${CMAKE_INSTALL_PREFIX}/include")
link_directories("${CMAKE_INSTALL_PREFIX}/lib")
