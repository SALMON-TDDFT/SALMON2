# CMake generated Testfile for 
# Source directory: /workspace/testsuites/303_classicEM_isolated_pulse_ttm
# Build directory: /workspace/build/testsuites/303_classicEM_isolated_pulse_ttm
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_303_classicEM_isolated_pulse_ttm "prep_303_classicEM_isolated_pulse_ttm.sh")
set_tests_properties(prep_303_classicEM_isolated_pulse_ttm PROPERTIES  FIXTURES_SETUP "setup_303_classicEM_isolated_pulse_ttm" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/303_classicEM_isolated_pulse_ttm/CMakeLists.txt;1;create_test;/workspace/testsuites/303_classicEM_isolated_pulse_ttm/CMakeLists.txt;0;")
add_test(run_303_classicEM_isolated_pulse_ttm "mpiexec_303_classicEM_isolated_pulse_ttm.sh")
set_tests_properties(run_303_classicEM_isolated_pulse_ttm PROPERTIES  FIXTURES_REQUIRED "setup_303_classicEM_isolated_pulse_ttm" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/303_classicEM_isolated_pulse_ttm/CMakeLists.txt;1;create_test;/workspace/testsuites/303_classicEM_isolated_pulse_ttm/CMakeLists.txt;0;")
add_test(verify_303_classicEM_isolated_pulse_ttm "verification")
set_tests_properties(verify_303_classicEM_isolated_pulse_ttm PROPERTIES  FIXTURES_CLEANUP "cleanup_303_classicEM_isolated_pulse_ttm" FIXTURES_REQUIRED "setup_303_classicEM_isolated_pulse_ttm;run_303_classicEM_isolated_pulse_ttm" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/303_classicEM_isolated_pulse_ttm/CMakeLists.txt;1;create_test;/workspace/testsuites/303_classicEM_isolated_pulse_ttm/CMakeLists.txt;0;")
