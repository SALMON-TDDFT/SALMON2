include(ExternalProject)

if (NOT USE_SCALAPACK)
  message(FATAL_ERROR "EigenExa requires ScaLAPACK library.")
endif ()

set(EIGENEXA_VERSION "2.4b")
message(STATUS "Build RIKEN R-CCS EigenExa library version ${EIGENEXA_VERSION}")

# set dummy option to avoid autotools checking...
set(HOST_OPTION)
if (${CMAKE_CROSSCOMPILING})
  set(HOST_OPTION "--host=x86_64")
endif ()

# compiler option, but skip fujitsu environment.
if (IS_FUJITSU_COMPILER)
else ()
  set(CC     "CC=${CMAKE_C_COMPILER}")
  set(FC     "F77=${CMAKE_Fortran_COMPILER}")
  set(CFLAGS "CFLAGS=${CMAKE_C_FLAGS}")
  set(FFLAGS "FFLAGS=${CMAKE_Fortran_FLAGS}")
endif()

include(${CMAKE_SOURCE_DIR}/cmakefiles/Builder/build_scalapack.cmake)

ExternalProject_Add(eigenexa-project
  URL               "https://www.r-ccs.riken.jp/labs/lpnctrt/assets/img/EigenExa-${EIGENEXA_VERSION}.tgz"
  PREFIX            ${CMAKE_BINARY_DIR}/eigenexa
  CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/eigenexa/src/eigenexa-project/configure ${HOST_OPTION} --prefix=${CMAKE_CURRENT_BINARY_DIR} ${CC} ${FC} ${CFLAGS} ${FFLAGS}
  BUILD_COMMAND     gmake -C src
  STEP_TARGETS      build
  EXCLUDE_FROM_ALL  on
)

add_library(eigenexa STATIC IMPORTED)
set_target_properties(eigenexa PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/eigenexa/src/eigenexa-project-build/src/libEigenExa.a)
add_dependencies(eigenexa eigenexa-project-build)
set(EXTERNAL_LIBS eigenexa ${EXTERNAL_LIBS})

# NOTE: EigenExa uses Fortran module files in src directory.
#       The compiler need that they can find the modules...
include_directories(${CMAKE_CURRENT_BINARY_DIR}/eigenexa/src/eigenexa-project-build/src)
