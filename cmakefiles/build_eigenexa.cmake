include(ExternalProject)

if (USE_EIGENEXA)
  if (NOT USE_SCALAPACK)
    message(FATAL_ERROR "EigenExa requires LAPACK and ScaLAPACK library.")
  endif ()

  set(EIGENEXA_VERSION "2.4b")
  message(STATUS "Build RIKEN R-CCS EigenExa library version ${EIGENEXA_VERSION}")

  ExternalProject_Add(eigenexa-project
    URL               "https://www.r-ccs.riken.jp/labs/lpnctrt/assets/img/EigenExa-${EIGENEXA_VERSION}.tgz"
    SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/eigenexa
    BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/eigenexa
    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/eigenexa/configure --prefix=${CMAKE_CURRENT_BINARY_DIR} CC=${CMAKE_C_COMPILER} F77=${CMAKE_Fortran_COMPILER} "CFLAGS=${CMAKE_C_FLAGS}" "FFLAGS=${CMAKE_Fortran_FLAGS}"
    BUILD_COMMAND     make
    STEP_TARGETS      install
    EXCLUDE_FROM_ALL  on
  )

  add_library(eigenexa STATIC IMPORTED)
  set_target_properties(eigenexa PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/eigenexa/src/libEigenExa.a)
  add_dependencies(eigenexa eigenexa-project-install)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} eigenexa)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/eigenexa/src)
endif ()
