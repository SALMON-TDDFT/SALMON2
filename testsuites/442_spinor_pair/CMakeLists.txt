add_executable(spinor_pair_probe_442
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_backend.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_exchange.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_ace.f90
  ${PROJECT_SOURCE_DIR}/samples/hse_mlwf_reference/native_spinor_pair_probe.f90)
set_target_properties(spinor_pair_probe_442 PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes/modules-442")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes/modules-442")
target_include_directories(spinor_pair_probe_442 PRIVATE "${CMAKE_BINARY_DIR}/testsuite-probes/modules-442")
target_link_libraries(spinor_pair_probe_442 PRIVATE ${EXTERNAL_LIBS})
if(EXTERNAL_PROJECT_TARGETS)
  add_dependencies(spinor_pair_probe_442 ${EXTERNAL_PROJECT_TARGETS})
endif()
set(TEST_COMMAND "${CMAKE_BINARY_DIR}/testsuite-probes/spinor_pair_probe_442")
create_test()
set_tests_properties(run_442_spinor_pair PROPERTIES FIXTURES_SETUP run_442_spinor_pair)
set_tests_properties(run_442_spinor_pair PROPERTIES TIMEOUT 120 PROCESSORS 1)
set(_pair_tests prep_442_spinor_pair run_442_spinor_pair verify_442_spinor_pair)
add_test(NAME 442_spinor_pair_large COMMAND $<TARGET_FILE:spinor_pair_probe_442> large)
set_tests_properties(442_spinor_pair_large PROPERTIES TIMEOUT 120 PROCESSORS 1)
list(APPEND _pair_tests 442_spinor_pair_large)
set_tests_properties(${_pair_tests} PROPERTIES
  LABELS "hybrid;spinor-pair"
  ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
