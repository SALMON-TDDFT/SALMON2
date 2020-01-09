### Numerical library
if (NOT ${CMAKE_CROSSCOMPILING})
  if (USE_SCALAPACK)
    find_package(ScaLAPACK REQUIRED)
    set(ScaLAPACK_FLAGS ${ScaLAPACK_LINKER_FLAGS})
  else ()
    find_package(LAPACK REQUIRED)
    set(LAPACK_FLAGS ${LAPACK_LINKER_FLAGS})
  endif ()
endif ()

if (USE_SCALAPACK)
  set(FILE_MATHLIB scalapack)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${ScaLAPACK_LIBRARIES})
  set(EXTERNAL_FLAGS ${EXTERNAL_FLAGS} ${ScaLAPACK_FLAGS})
else ()
  set(FILE_MATHLIB lapack)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${LAPACK_LIBRARIES})
  set(EXTERNAL_FLAGS ${EXTERNAL_FLAGS} ${LAPACK_FLAGS})
endif ()


### Libxc
if (USE_LIBXC)
  find_package(Libxc REQUIRED)
  include_directories($<TARGET_PROPERTY:Libxc::xc,INTERFACE_INCLUDE_DIRECTORIES>)
  link_libraries     (Libxc::xcf90 Libxc::xc)
endif ()
