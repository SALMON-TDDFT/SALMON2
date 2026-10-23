add_executable(spinor_pair_probe_443
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_backend.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_collinear_fft.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_exchange.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor.f90
  ${PROJECT_SOURCE_DIR}/samples/hse_mlwf_reference/native_spinor_pair_mpi_probe.f90)
set_target_properties(spinor_pair_probe_443 PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes/modules-443")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes/modules-443")
target_include_directories(spinor_pair_probe_443 BEFORE PRIVATE "${CMAKE_BINARY_DIR}/testsuite-probes/modules-443")
target_link_libraries(spinor_pair_probe_443 PRIVATE ${EXTERNAL_LIBS})
if(EXTERNAL_PROJECT_TARGETS)
  add_dependencies(spinor_pair_probe_443 ${EXTERNAL_PROJECT_TARGETS})
endif()
set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_BINARY_DIR}/testsuite-probes/spinor_pair_probe_443")
create_mpi_test()
set_tests_properties(run_443_spinor_pair_mpi PROPERTIES FIXTURES_SETUP run_443_spinor_pair_mpi)
set_tests_properties(run_443_spinor_pair_mpi PROPERTIES TIMEOUT 120 PROCESSORS 4)
set(_pair_tests prep_443_spinor_pair_mpi run_443_spinor_pair_mpi verify_443_spinor_pair_mpi)
add_test(NAME 443_spinor_pair_mpi_empty COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 $<TARGET_FILE:spinor_pair_probe_443>)
set_tests_properties(443_spinor_pair_mpi_empty PROPERTIES TIMEOUT 120 PROCESSORS 8)
list(APPEND _pair_tests 443_spinor_pair_mpi_empty)
add_test(NAME 443_spinor_pair_mpi_contiguous COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:spinor_pair_probe_443> contiguous)
set_tests_properties(443_spinor_pair_mpi_contiguous PROPERTIES TIMEOUT 120 PROCESSORS 4)
list(APPEND _pair_tests 443_spinor_pair_mpi_contiguous)
add_test(NAME 443_spinor_pair_mpi_nan COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:spinor_pair_probe_443> nan)
set_tests_properties(443_spinor_pair_mpi_nan PROPERTIES TIMEOUT 120 PROCESSORS 4)
list(APPEND _pair_tests 443_spinor_pair_mpi_nan)
add_test(NAME 443_spinor_pair_mpi_mixed COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:spinor_pair_probe_443> mixed)
set_tests_properties(443_spinor_pair_mpi_mixed PROPERTIES TIMEOUT 120 PROCESSORS 4)
list(APPEND _pair_tests 443_spinor_pair_mpi_mixed)
add_test(NAME 443_spinor_pair_mpi_overflow COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:spinor_pair_probe_443> overflow)
set_tests_properties(443_spinor_pair_mpi_overflow PROPERTIES TIMEOUT 120 PROCESSORS 4)
list(APPEND _pair_tests 443_spinor_pair_mpi_overflow)
set_tests_properties(${_pair_tests} PROPERTIES
  LABELS "hybrid;spinor-pair"
  ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
