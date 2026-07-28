include(ExternalProject)
find_program(SALMON_EXTERNAL_MAKE NAMES gmake make REQUIRED)

if ((CMAKE_CROSSCOMPILING OR IS_FUJITSU_COMPILER) AND
    NOT SALMON_PETSC_CROSS_BUILD)
  message(FATAL_ERROR
    "Automatic PETSc/SLEPc builds are not supported for this cross or Fujitsu build. "
    "Install PETSc and SLEPc for the target platform and set SLEPC_INSTALLDIR "
    "(and PETSC_INSTALLDIR when needed).")
endif ()

set(PETSC_VERSION "3.25.3")
set(SALMON_PETSC_PREFIX "${CMAKE_BINARY_DIR}/petsc-install")

if (SALMON_PETSC_COPTFLAGS)
  set(_petsc_coptflags "${SALMON_PETSC_COPTFLAGS}")
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(_petsc_coptflags "${CMAKE_C_FLAGS_DEBUG}")
else ()
  set(_petsc_coptflags "${CMAKE_C_FLAGS_RELEASE}")
endif ()

if (SALMON_PETSC_FOPTFLAGS)
  set(_petsc_foptflags "${SALMON_PETSC_FOPTFLAGS}")
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(_petsc_foptflags "${CMAKE_Fortran_FLAGS_DEBUG}")
else ()
  set(_petsc_foptflags "${CMAKE_Fortran_FLAGS_RELEASE}")
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(_petsc_debugging 1)
else ()
  set(_petsc_debugging 0)
endif ()

set(_petsc_blaslapack_options)
if (LAPACK_INSTALLDIR)
  list(APPEND _petsc_blaslapack_options
       "--with-blaslapack-dir=${LAPACK_INSTALLDIR}")
elseif (SALMON_BLAS_LIBRARY AND SALMON_LAPACK_LIBRARY)
  list(APPEND _petsc_blaslapack_options
       "--with-blas-lib=${SALMON_BLAS_LIBRARY}"
       "--with-lapack-lib=${SALMON_LAPACK_LIBRARY}")
elseif (SALMON_BLASLAPACK_LIBRARIES)
  set(_petsc_uses_accelerate FALSE)
  if (APPLE)
    foreach (_blaslapack_library IN LISTS SALMON_BLASLAPACK_LIBRARIES)
      if (_blaslapack_library MATCHES "Accelerate\\.framework")
        set(_petsc_uses_accelerate TRUE)
      endif ()
    endforeach ()
  endif ()
  if (_petsc_uses_accelerate)
    list(APPEND _petsc_blaslapack_options
         "--with-blaslapack-lib=[-framework,Accelerate]")
  else ()
    list(JOIN SALMON_BLASLAPACK_LIBRARIES "," _petsc_blaslapack_libraries)
    list(APPEND _petsc_blaslapack_options
         "--with-blaslapack-lib=[${_petsc_blaslapack_libraries}]")
  endif ()
else ()
  message(FATAL_ERROR
    "The BLAS/LAPACK libraries selected for SALMON could not be passed to PETSc.")
endif ()

set(_petsc_extra_options)
if (PETSC_CONFIGURE_OPTIONS)
  separate_arguments(_petsc_extra_options NATIVE_COMMAND
                     "${PETSC_CONFIGURE_OPTIONS}")
endif ()

set(_petsc_platform_options)
if (SALMON_PETSC_PLATFORM_OPTIONS)
  separate_arguments(_petsc_platform_options NATIVE_COMMAND
                     "${SALMON_PETSC_PLATFORM_OPTIONS}")
endif ()

set(_petsc_dependencies)
if (SALMON_LAPACK_BUILD_TARGET)
  list(APPEND _petsc_dependencies ${SALMON_LAPACK_BUILD_TARGET})
endif ()

set(_petsc_patch_command)
if (IS_FUJITSU_COMPILER)
  set(_petsc_patch_command
      ${CMAKE_COMMAND}
      "-DPETSC_SOURCE_DIR=<SOURCE_DIR>"
      -P "${CMAKE_SOURCE_DIR}/cmakefiles/Builder/patch_petsc_fujitsu.cmake")
endif ()

message(STATUS "Build PETSc version ${PETSC_VERSION}")
ExternalProject_Add(petsc-project
  URL               "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PETSC_VERSION}.tar.gz"
  URL_HASH          SHA256=95ce60df2c7f9c5044d6a544c41e996a512557f91df1a60bdb690b332904ebb5
  PREFIX            "${CMAKE_BINARY_DIR}/petsc"
  BUILD_IN_SOURCE   on
  PATCH_COMMAND     ${_petsc_patch_command}
  CONFIGURE_COMMAND <SOURCE_DIR>/configure
                    "--prefix=${SALMON_PETSC_PREFIX}"
                    "--with-cc=${CMAKE_C_COMPILER}"
                    "--with-cxx=0"
                    "--with-fc=${CMAKE_Fortran_COMPILER}"
                    "--with-mpi=1"
                    "--with-x=0"
                    "--with-debugging=${_petsc_debugging}"
                    "--with-shared-libraries=0"
                    "--with-scalar-type=real"
                    "--with-precision=double"
                    "--with-64-bit-indices=0"
                    "COPTFLAGS=${_petsc_coptflags}"
                    "FOPTFLAGS=${_petsc_foptflags}"
                    ${_petsc_blaslapack_options}
                    ${_petsc_platform_options}
                    ${_petsc_extra_options}
  BUILD_COMMAND     ${CMAKE_COMMAND} -E env
                    "PETSC_DIR=<SOURCE_DIR>" ${SALMON_EXTERNAL_MAKE} all
  INSTALL_COMMAND   ${CMAKE_COMMAND} -E env
                    "PETSC_DIR=<SOURCE_DIR>" ${SALMON_EXTERNAL_MAKE} install
  DEPENDS           ${_petsc_dependencies}
  STEP_TARGETS      install
  EXCLUDE_FROM_ALL  on
)

add_library(petsc STATIC IMPORTED GLOBAL)
set_target_properties(petsc PROPERTIES
  IMPORTED_LOCATION "${SALMON_PETSC_PREFIX}/lib/libpetsc.a")
add_dependencies(petsc petsc-project-install)
