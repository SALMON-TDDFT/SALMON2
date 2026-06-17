include(ExternalProject)

set(WANNIER90_VERSION "3.1.0" CACHE STRING "Wannier90 version")
message(STATUS "Build Wannier90 library version ${WANNIER90_VERSION}")

find_program(WANNIER90_MAKE_EXECUTABLE NAMES gmake make)
if (NOT WANNIER90_MAKE_EXECUTABLE)
  message(FATAL_ERROR "Wannier90 build requires make or gmake.")
endif ()

set(WANNIER90_PREFIX ${CMAKE_BINARY_DIR}/wannier90)
set(WANNIER90_MAKE_INC ${WANNIER90_PREFIX}/make.inc)
set(WANNIER90_LIBPATH ${WANNIER90_PREFIX}/src/wannier90-project/libwannier.a)
set(WANNIER90_EXECUTABLE_PATH ${WANNIER90_PREFIX}/src/wannier90-project/wannier90.x CACHE FILEPATH "Wannier90 executable")
set(WANNIER90_COMMS "serial")
set(WANNIER90_BUILD_TARGETS lib)
if (USE_MPI)
  set(WANNIER90_COMMS "mpi")
  set(WANNIER90_BUILD_TARGETS wannier lib)
endif ()
set(WANNIER90_COMMS "${WANNIER90_COMMS}" CACHE STRING "Wannier90 communication backend")
set_property(CACHE WANNIER90_COMMS PROPERTY STRINGS serial mpi)
message(STATUS "Build Wannier90 with COMMS=${WANNIER90_COMMS}")
if (NOT LAPACK_LIBRARIES)
  find_package(LAPACK QUIET)
endif ()
set(WANNIER90_LINK_LIBS "${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}")
string(REPLACE ";" " " WANNIER90_LINK_LIBS "${WANNIER90_LINK_LIBS}")
if (APPLE AND WANNIER90_LINK_LIBS MATCHES "Accelerate\\.framework")
  set(WANNIER90_LINK_LIBS "-framework Accelerate -lm -ldl")
endif ()
if (APPLE AND NOT WANNIER90_LINK_LIBS)
  set(WANNIER90_LINK_LIBS "-framework Accelerate -lm -ldl")
endif ()

file(MAKE_DIRECTORY ${WANNIER90_PREFIX})
file(WRITE ${WANNIER90_MAKE_INC}
"F90 = ${CMAKE_Fortran_COMPILER}
COMMS = ${WANNIER90_COMMS}
MPIF90 = ${CMAKE_Fortran_COMPILER}
FCOPTS = ${CMAKE_Fortran_FLAGS}
LDOPTS = ${CMAKE_EXE_LINKER_FLAGS}
LIBS = ${WANNIER90_LINK_LIBS}

AR = ${CMAKE_AR}
RANLIB = ${CMAKE_RANLIB}

default: lib
")

ExternalProject_Add(wannier90-project
  URL               "https://github.com/wannier-developers/wannier90/archive/refs/tags/v${WANNIER90_VERSION}.tar.gz"
  PREFIX            ${WANNIER90_PREFIX}
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy ${WANNIER90_MAKE_INC} <SOURCE_DIR>/make.inc
  BUILD_COMMAND     ${WANNIER90_MAKE_EXECUTABLE} ${WANNIER90_BUILD_TARGETS}
  INSTALL_COMMAND   ""
  BUILD_IN_SOURCE   1
)

add_library(wannier90 STATIC IMPORTED)
set_target_properties(wannier90 PROPERTIES IMPORTED_LOCATION ${WANNIER90_LIBPATH})
add_dependencies(wannier90 wannier90-project)
set(EXTERNAL_LIBS wannier90 ${EXTERNAL_LIBS})
