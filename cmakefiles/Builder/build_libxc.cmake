# One Libxc serves both native HSE (C ABI) and optional conventional Libxc XC.
include(ExternalProject)
include(CheckCSourceCompiles)
include(CheckFortranSourceCompiles)
include(CMakePushCheckState)
include(${CMAKE_CURRENT_LIST_DIR}/hse_external_options.cmake)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(PC_SALMON_XC QUIET libxc)
endif()
find_path(SALMON_XC_INCLUDE xc.h HINTS ${LIBXC_INSTALLDIR}/include ${PC_SALMON_XC_INCLUDE_DIRS})
find_library(SALMON_XC_LIBRARY xc HINTS ${LIBXC_INSTALLDIR}/lib ${PC_SALMON_XC_LIBRARY_DIRS})
set(_xc_usable FALSE)
if(SALMON_XC_INCLUDE AND SALMON_XC_LIBRARY)
  cmake_push_check_state(RESET)
  set(CMAKE_REQUIRED_INCLUDES ${SALMON_XC_INCLUDE})
  set(CMAKE_REQUIRED_LIBRARIES ${SALMON_XC_LIBRARY} m)
  unset(SALMON_XC_C_WORKS CACHE)
  check_c_source_compiles("#include <xc.h>
#include <xc_version.h>
#if XC_MAJOR_VERSION < 5
#error Libxc 5 or later required
#endif
int main(void) {
 xc_func_type *p = xc_func_alloc(); double w,a,b,v[3]={.25,.11,.11};
 if(xc_func_init(p,428,1)) return 1;
 xc_func_set_ext_params(p,v); xc_hyb_cam_coef(p,&w,&a,&b);
 xc_func_end(p); xc_func_free(p); return 0;
}" SALMON_XC_C_WORKS)
  set(_xc_usable ${SALMON_XC_C_WORKS})
  if(USE_LIBXC)
    find_library(SALMON_XCF90_LIBRARY xcf90 HINTS ${LIBXC_INSTALLDIR}/lib ${PC_SALMON_XC_LIBRARY_DIRS})
    find_path(SALMON_XCF90_INCLUDE xc_f90_lib_m.mod HINTS ${SALMON_XC_INCLUDE} ${LIBXC_INSTALLDIR}/include)
    set(_xc_usable FALSE)
    if(SALMON_XC_C_WORKS AND SALMON_XCF90_LIBRARY AND SALMON_XCF90_INCLUDE)
      set(CMAKE_REQUIRED_INCLUDES ${SALMON_XC_INCLUDE} ${SALMON_XCF90_INCLUDE})
      set(CMAKE_REQUIRED_LIBRARIES ${SALMON_XCF90_LIBRARY} ${SALMON_XC_LIBRARY} m)
      unset(SALMON_XC_FORTRAN_WORKS CACHE)
      check_fortran_source_compiles("program probe
use xc_f90_lib_m
implicit none
integer :: major,minor,micro
call xc_f90_version(major,minor,micro)
end program" SALMON_XC_FORTRAN_WORKS SRC_EXT F90)
      set(_xc_usable ${SALMON_XC_FORTRAN_WORKS})
    endif()
  endif()
  cmake_pop_check_state()
endif()
if(_xc_usable)
  message(STATUS "Using compatible installed Libxc: ${SALMON_XC_LIBRARY}")
  include_directories(${SALMON_XC_INCLUDE})
  if(USE_LIBXC)
    include_directories(${SALMON_XCF90_INCLUDE})
    list(APPEND EXTERNAL_LIBS ${SALMON_XCF90_LIBRARY})
  endif()
  list(APPEND EXTERNAL_LIBS ${SALMON_XC_LIBRARY} m)
else()
  # 5.2.3 provides both the C API and SALMON's legacy xc_f90 module interface.
  set(_prefix "${CMAKE_BINARY_DIR}/dependencies/libxc")
  message(STATUS "Automatically building compatible Libxc 5.2.3")
  set(_xc_byproducts "${_prefix}/lib/libxc.a")
  if(USE_LIBXC)
    list(APPEND _xc_byproducts "${_prefix}/lib/libxcf90.a" "${_prefix}/lib/libxcf03.a")
  endif()
  ExternalProject_Add(libxc-project
    URL https://gitlab.com/libxc/libxc/-/archive/5.2.3/libxc-5.2.3.tar.bz2
    URL_HASH SHA256=851a45aee9ddaafea49f684fe3e3ab2fbd79f1c1289b8c1cea216330b27e887b
    PREFIX "${CMAKE_BINARY_DIR}/dependencies/libxc-build"
    LIST_SEPARATOR "|"
    CMAKE_ARGS ${SALMON_EXTERNAL_CMAKE_ARGS} -DCMAKE_INSTALL_PREFIX:PATH=${_prefix}
      -DENABLE_FORTRAN:BOOL=${USE_LIBXC} -DENABLE_XHOST:BOOL=OFF -DBUILD_TESTING:BOOL=OFF
    BUILD_BYPRODUCTS ${_xc_byproducts})
  file(MAKE_DIRECTORY "${_prefix}/include")
  include_directories("${_prefix}/include")
  add_library(salmon_xc STATIC IMPORTED)
  set_target_properties(salmon_xc PROPERTIES IMPORTED_LOCATION "${_prefix}/lib/libxc.a")
  add_dependencies(salmon_xc libxc-project)
  if(USE_LIBXC)
    add_library(salmon_xcf90 STATIC IMPORTED)
    set_target_properties(salmon_xcf90 PROPERTIES IMPORTED_LOCATION "${_prefix}/lib/libxcf90.a")
    add_dependencies(salmon_xcf90 libxc-project)
    list(APPEND EXTERNAL_LIBS salmon_xcf90)
  endif()
  list(APPEND EXTERNAL_LIBS salmon_xc m)
  list(APPEND EXTERNAL_PROJECT_TARGETS libxc-project)
endif()
