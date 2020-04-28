include(ExternalProject)

find_package(Libxc QUIET)

if (Libxc_FOUND)
  message(STATUS "Libxc found.")
  include_directories($<TARGET_PROPERTY:Libxc::xc,INTERFACE_INCLUDE_DIRECTORIES>)
  link_libraries(Libxc::xcf90 Libxc::xc)
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
  set(EXTERNAL_LIBS xcf90 xc ${EXTERNAL_LIBS})

  # NOTE: Libxc will installs the Fortran module files to `include` directory.
  #       The compiler need that they can find the modules...
  include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
endif ()
