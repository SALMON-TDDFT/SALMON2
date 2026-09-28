if (LAPACK_INSTALLDIR)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${LAPACK_INSTALLDIR})
endif ()

if (USE_EIGENEXA)
  # Find or Build EigenExa, ScaLAPACK and LAPACK
  include(${CMAKE_SOURCE_DIR}/cmakefiles/Builder/build_eigenexa.cmake)
elseif (USE_SCALAPACK)
  # Find or Build ScaLAPACK and LAPACK
  include(${CMAKE_SOURCE_DIR}/cmakefiles/Builder/build_scalapack.cmake)
else ()
  # Find or Build LAPACK
  include(${CMAKE_SOURCE_DIR}/cmakefiles/Builder/build_lapack.cmake)
endif ()

if (USE_LIBXC)
  if (LIBXC_INSTALLDIR)
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${LIBXC_INSTALLDIR})
  endif ()
  include(${CMAKE_SOURCE_DIR}/cmakefiles/Builder/build_libxc.cmake)
endif ()
