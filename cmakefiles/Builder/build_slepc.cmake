include(ExternalProject)

if (SLEPC_INSTALLDIR)
  list(APPEND CMAKE_PREFIX_PATH ${SLEPC_INSTALLDIR})
endif ()
if (PETSC_INSTALLDIR)
  list(APPEND CMAKE_PREFIX_PATH ${PETSC_INSTALLDIR})
endif ()

find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules(SLEPC QUIET IMPORTED_TARGET GLOBAL SLEPc)
endif ()

if (SLEPC_FOUND)
  message(STATUS "Found SLEPc version ${SLEPC_VERSION}")
  set(EXTERNAL_LIBS PkgConfig::SLEPC ${EXTERNAL_LIBS})
else ()
  set(_slepc_search_paths ${SLEPC_INSTALLDIR})
  if (DEFINED ENV{SLEPC_DIR})
    list(APPEND _slepc_search_paths "$ENV{SLEPC_DIR}")
  endif ()
  set(_petsc_search_paths ${PETSC_INSTALLDIR})
  if (DEFINED ENV{PETSC_DIR})
    list(APPEND _petsc_search_paths "$ENV{PETSC_DIR}")
  endif ()

  find_path(SLEPC_INCLUDE_DIR slepceps.mod
            PATHS ${_slepc_search_paths} PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_library(SLEPC_LIBRARY NAMES slepc
               PATHS ${_slepc_search_paths} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
  find_path(PETSC_INCLUDE_DIR petscsys.mod
            PATHS ${_petsc_search_paths} ${_slepc_search_paths}
            PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_library(PETSC_LIBRARY NAMES petsc
               PATHS ${_petsc_search_paths} ${_slepc_search_paths}
               PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)

  if (SLEPC_INCLUDE_DIR AND SLEPC_LIBRARY AND
      PETSC_INCLUDE_DIR AND PETSC_LIBRARY)
    message(STATUS "Found SLEPc library: ${SLEPC_LIBRARY}")
    include_directories(${SLEPC_INCLUDE_DIR} ${PETSC_INCLUDE_DIR})
    set(EXTERNAL_LIBS ${SLEPC_LIBRARY} ${PETSC_LIBRARY} ${EXTERNAL_LIBS})
  else ()
    include(${CMAKE_SOURCE_DIR}/cmakefiles/Builder/build_petsc.cmake)

    set(SLEPC_VERSION "3.25.1")
    set(SALMON_SLEPC_PREFIX "${CMAKE_BINARY_DIR}/slepc-install")
    message(STATUS "Build SLEPc version ${SLEPC_VERSION}")

    ExternalProject_Add(slepc-project
      URL               "https://slepc.upv.es/download/distrib/slepc-${SLEPC_VERSION}.tar.gz"
      URL_HASH          SHA256=906ddbe15a20774c23ddcdf13a5054889d00a26c3c37463447ee593c757d03ee
      PREFIX            "${CMAKE_BINARY_DIR}/slepc"
      BUILD_IN_SOURCE   on
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env
                        "PETSC_DIR=${SALMON_PETSC_PREFIX}"
                        <SOURCE_DIR>/configure
                        "--prefix=${SALMON_SLEPC_PREFIX}"
      BUILD_COMMAND     ${CMAKE_COMMAND} -E env
                        "PETSC_DIR=${SALMON_PETSC_PREFIX}"
                        "SLEPC_DIR=<SOURCE_DIR>" ${SALMON_EXTERNAL_MAKE}
      INSTALL_COMMAND   ${CMAKE_COMMAND} -E env
                        "PETSC_DIR=${SALMON_PETSC_PREFIX}"
                        "SLEPC_DIR=<SOURCE_DIR>" ${SALMON_EXTERNAL_MAKE} install
      DEPENDS           petsc-project-install
      STEP_TARGETS      install
      EXCLUDE_FROM_ALL  on
    )

    add_library(slepc STATIC IMPORTED GLOBAL)
    set_target_properties(slepc PROPERTIES
      IMPORTED_LOCATION "${SALMON_SLEPC_PREFIX}/lib/libslepc.a")
    add_dependencies(slepc slepc-project-install)

    include_directories("${SALMON_SLEPC_PREFIX}/include"
                        "${SALMON_PETSC_PREFIX}/include")
    find_package(Threads REQUIRED)
    set(EXTERNAL_LIBS slepc petsc ${EXTERNAL_LIBS}
                      Threads::Threads ${CMAKE_DL_LIBS})
  endif ()
endif ()
