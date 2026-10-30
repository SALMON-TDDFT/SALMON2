add_executable(collinear_worker_probe
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_backend.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_collinear_fft.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_exchange.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_symmetry.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor_wave_fft.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_worker_layout.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor_wave_workers.f90 worker_probe.f90)
set(_worker_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-462-worker")
file(MAKE_DIRECTORY "${_worker_modules}")
set_target_properties(collinear_worker_probe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${_worker_modules}")
target_include_directories(collinear_worker_probe BEFORE PRIVATE "${_worker_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(collinear_worker_probe PRIVATE ${EXTERNAL_LIBS})
salmon_guard_eigensolver_target(collinear_worker_probe)
foreach(_np 2 4 8 16 33)
 add_test(NAME 462_collinear_worker_${_np} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${_np}
  $<TARGET_FILE:collinear_worker_probe>)
 set_tests_properties(462_collinear_worker_${_np} PROPERTIES TIMEOUT 60 PROCESSORS ${_np}
  LABELS "hybrid;collinear;openmp;mpi" ENVIRONMENT "OMP_NUM_THREADS=4;OMP_DYNAMIC=FALSE;OPENBLAS_NUM_THREADS=1")
endforeach()

add_executable(worker_layout_probe ${PROJECT_SOURCE_DIR}/src/xc/exx_worker_layout.f90 layout_probe.f90)
set(_layout_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-462-layout")
file(MAKE_DIRECTORY "${_layout_modules}")
set_target_properties(worker_layout_probe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${_layout_modules}")
salmon_guard_eigensolver_target(worker_layout_probe)

add_executable(balanced_strip_probe ${PROJECT_SOURCE_DIR}/src/xc/exx_worker_layout.f90 balanced_probe.f90)
set(_balanced_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-462-balanced")
file(MAKE_DIRECTORY "${_balanced_modules}")
set_target_properties(balanced_strip_probe PROPERTIES
 RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes" Fortran_MODULE_DIRECTORY "${_balanced_modules}")
salmon_guard_eigensolver_target(balanced_strip_probe)

# Python-based development checks are opt-in and absent from public source archives.
option(SALMON_ENABLE_DEVELOPMENT_TESTS "Enable development regression checks" OFF)
set(_development_tests "${PROJECT_SOURCE_DIR}/tools/development-tests/462_collinear_workers")
if(SALMON_ENABLE_DEVELOPMENT_TESTS)
 if(NOT EXISTS "${_development_tests}/CMakeLists.txt")
  message(FATAL_ERROR "Development tests are unavailable in this source distribution")
 endif()
 add_subdirectory("${_development_tests}" "${CMAKE_CURRENT_BINARY_DIR}/development")
endif()
