create_test()
set_tests_properties(run_142_Si_tdcdft PROPERTIES
  FIXTURES_REQUIRED "setup_142_Si_tdcdft;run_111_bulk_Si_gs_dp")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/unit_modules)
add_executable(test_tdcdft_field ${CMAKE_SOURCE_DIR}/src/rt/tdcdft_lrc.f90 test_update.f90)
set_target_properties(test_tdcdft_field PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test_bins
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/unit_modules)
target_include_directories(test_tdcdft_field PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/unit_modules)
add_test(NAME unit_tdcdft_field COMMAND test_tdcdft_field)
