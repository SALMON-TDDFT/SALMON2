cmake_minimum_required(VERSION 3.14)
project(SAWFSpglibBackendTest LANGUAGES C Fortran)

if (NOT SALMON_ROOT)
  message(FATAL_ERROR "SALMON_ROOT is required")
endif ()

option(SAWF_TEST_HAVE_SPGLIB "Build the focused SAWF test with Spglib" OFF)
set(HAVE_SPGLIB ${SAWF_TEST_HAVE_SPGLIB})
configure_file(config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
configure_file(
  ${SALMON_ROOT}/src/gs/dc/lcfo_wannier_sawf.f90
  ${CMAKE_CURRENT_BINARY_DIR}/lcfo_wannier_sawf.F90
  COPYONLY
)

set(SAWF_TEST_SOURCES
  mock_sym_sub.f90
  ${CMAKE_CURRENT_BINARY_DIR}/lcfo_wannier_sawf.F90
  driver.F90
)
if (SAWF_TEST_HAVE_SPGLIB)
  find_package(Spglib 2.1 CONFIG REQUIRED)
  list(APPEND SAWF_TEST_SOURCES
    ${SALMON_ROOT}/src/gs/dc/lcfo_wannier_spglib.c
  )
endif ()

add_executable(check_sawf_spglib_backend ${SAWF_TEST_SOURCES})
target_include_directories(check_sawf_spglib_backend PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
)
if (SAWF_TEST_HAVE_SPGLIB)
  target_link_libraries(check_sawf_spglib_backend PRIVATE Spglib::symspg)

  add_executable(check_sawf_spglib_capacity
    ${SALMON_ROOT}/src/gs/dc/lcfo_wannier_spglib.c
    check_wrapper_capacity.c
  )
  target_link_libraries(check_sawf_spglib_capacity PRIVATE Spglib::symspg)
endif ()

enable_testing()
add_test(NAME sawf_spglib_backend COMMAND check_sawf_spglib_backend)
if (SAWF_TEST_HAVE_SPGLIB)
  add_test(NAME sawf_spglib_capacity COMMAND check_sawf_spglib_capacity)
endif ()
