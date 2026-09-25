include(ExternalProject)
include(${CMAKE_CURRENT_LIST_DIR}/hse_external_options.cmake)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(PC_HSE_FFTW QUIET fftw3)
endif()
find_path(HSE_FFTW_INCLUDE fftw3.f03 HINTS ${FFTW_INSTALLDIR}/include ${PC_HSE_FFTW_INCLUDE_DIRS})
find_library(HSE_FFTW_LIBRARY fftw3 HINTS ${FFTW_INSTALLDIR}/lib ${PC_HSE_FFTW_LIBRARY_DIRS})
include(CheckFortranSourceCompiles)
include(CMakePushCheckState)
set(_fftw_usable FALSE)
if(HSE_FFTW_INCLUDE AND HSE_FFTW_LIBRARY)
  cmake_push_check_state(RESET)
  set(CMAKE_REQUIRED_INCLUDES ${HSE_FFTW_INCLUDE})
  set(CMAKE_REQUIRED_LIBRARIES ${HSE_FFTW_LIBRARY} m)
  unset(SALMON_FFTW_WORKS CACHE)
  check_fortran_source_compiles("program probe
use iso_c_binding
implicit none
include 'fftw3.f03'
type(c_ptr) :: p
p=fftw_alloc_complex(1_c_size_t)
call fftw_free(p)
end program" SALMON_FFTW_WORKS SRC_EXT F90)
  set(_fftw_usable ${SALMON_FFTW_WORKS})
  cmake_pop_check_state()
endif()
if(_fftw_usable)
  message(STATUS "HSE: using installed FFTW ${HSE_FFTW_LIBRARY}")
  include_directories(${HSE_FFTW_INCLUDE})
  list(APPEND EXTERNAL_LIBS ${HSE_FFTW_LIBRARY})
else()
  set(_prefix "${CMAKE_BINARY_DIR}/dependencies/fftw")
  message(STATUS "HSE: automatically building FFTW 3.3.10")
  ExternalProject_Add(hse-fftw-project
    URL https://www.fftw.org/fftw-3.3.10.tar.gz
    URL_HASH SHA256=56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467
    PREFIX "${CMAKE_BINARY_DIR}/dependencies/fftw-build"
    LIST_SEPARATOR "|"
    CMAKE_ARGS ${SALMON_EXTERNAL_CMAKE_ARGS} -DCMAKE_INSTALL_PREFIX:PATH=${_prefix}
      -DBUILD_TESTS:BOOL=OFF -DENABLE_OPENMP:BOOL=OFF -DENABLE_THREADS:BOOL=OFF
    BUILD_BYPRODUCTS "${_prefix}/lib/libfftw3.a")
  file(MAKE_DIRECTORY "${_prefix}/include")
  add_library(salmon_hse_fftw STATIC IMPORTED)
  set_target_properties(salmon_hse_fftw PROPERTIES IMPORTED_LOCATION "${_prefix}/lib/libfftw3.a")
  add_dependencies(salmon_hse_fftw hse-fftw-project)
  list(APPEND EXTERNAL_PROJECT_TARGETS hse-fftw-project)
  include_directories("${_prefix}/include")
  list(APPEND EXTERNAL_LIBS salmon_hse_fftw)
endif()
