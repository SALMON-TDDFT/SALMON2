add_executable(exx_k_selection_probe
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_selection.f90 ${CMAKE_CURRENT_SOURCE_DIR}/probe.f90)
set_target_properties(exx_k_selection_probe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes/modules-445")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes/modules-445")
target_include_directories(exx_k_selection_probe PRIVATE "${CMAKE_BINARY_DIR}/testsuite-probes/modules-445")
set(TEST_COMMAND "${CMAKE_BINARY_DIR}/testsuite-probes/exx_k_selection_probe")
create_test()
set_tests_properties(run_445_exx_k_selection PROPERTIES FIXTURES_SETUP run_445_exx_k_selection TIMEOUT 30)
set_tests_properties(prep_445_exx_k_selection run_445_exx_k_selection verify_445_exx_k_selection
  PROPERTIES LABELS "hybrid;exchange-selection")
add_test(NAME 445_wse2_preparation COMMAND python3 -m unittest discover
  -s ${PROJECT_SOURCE_DIR}/samples/wse2_hse_soc_nlcc -p test_preparation.py)
set_tests_properties(445_wse2_preparation PROPERTIES LABELS "hybrid;exchange-selection" TIMEOUT 30)
