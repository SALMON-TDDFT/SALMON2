include(ExternalProject)


### LAPACK
find_package(LAPACK QUIET)

if (LAPACK_FOUND)
  message(STATUS "LAPACK library found.")
else ()
  # NOTE: LAPACK 3.7.0 and later version can't build by GCC 4.8.5, which RHEL7 provided compiler.
  set(LAPACK_VERSION "3.6.1")
  message(STATUS "Build Netlib LAPACK library version ${LAPACK_VERSION}")

  ExternalProject_Add(lapack-project
    URL              "http://www.netlib.org/lapack/lapack-${LAPACK_VERSION}.tgz"
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
  set(EXTERNAL_LIBS lapack blas ${EXTERNAL_LIBS})
endif ()


### ScaLAPACK
if (USE_SCALAPACK)
  if (USE_MPI)
  else ()
    message(FATAL_ERROR "Use ScaLAPACK: but MPI feature disabled.")
  endif ()

  find_package(ScaLAPACK QUIET)

  if (ScaLAPACK_FOUND)
    message(STATUS "ScaLAPACK library found.")
  else ()
    set(SCALAPACK_VERSION "2.1.0")
    message(STATUS "Build Netlib ScaLAPACK library version ${SCALAPACK_VERSION}")

    ExternalProject_Add(scalapack-project
      URL              "https://github.com/Reference-ScaLAPACK/scalapack/archive/v${SCALAPACK_VERSION}.tar.gz"
      PREFIX           "${CMAKE_BINARY_DIR}/scalapack"
      CMAKE_ARGS       -D BUILD_SHARED_LIBS=off
                       -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
                       -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER} -D CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                       -D CMAKE_C_FLAGS=${CMAKE_C_FLAGS} -D CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
                       -D CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG} -D CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}
                       -D CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE} -D CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}
      STEP_TARGETS     install
      EXCLUDE_FROM_ALL on
    )

    add_library(scalapack STATIC IMPORTED)
    set_target_properties(scalapack PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib/libscalapack.a)
    add_dependencies(scalapack-project-install    lapack-project-install)
    add_dependencies(scalapack                 scalapack-project-install)
    set(EXTERNAL_LIBS scalapack ${EXTERNAL_LIBS})
  endif ()
endif ()


if (USE_SCALAPACK)
  if (ScaLAPACK_FOUND)
    set(EXTERNAL_LIBS  ${EXTERNAL_LIBS}  ${ScaLAPACK_LIBRARIES})
    set(EXTERNAL_FLAGS ${EXTERNAL_FLAGS} ${ScaLAPACK_FLAGS})
  endif ()
else()
  if (LAPACK_FOUND)
    set(EXTERNAL_LIBS  ${EXTERNAL_LIBS}  ${LAPACK_LIBRARIES})
    set(EXTERNAL_FLAGS ${EXTERNAL_FLAGS} ${LAPACK_FLAGS})
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
      URL              "http://www.tddft.org/programs/libxc/down.php?file=${LIBXC_VERSION}/libxc-${LIBXC_VERSION}.tar.gz"
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
