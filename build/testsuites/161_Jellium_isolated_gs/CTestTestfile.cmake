# CMake generated Testfile for 
# Source directory: /workspace/testsuites/161_Jellium_isolated_gs
# Build directory: /workspace/build/testsuites/161_Jellium_isolated_gs
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_161_Jellium_isolated_gs "prep_161_Jellium_isolated_gs.sh")
set_tests_properties(prep_161_Jellium_isolated_gs PROPERTIES  FIXTURES_SETUP "setup_161_Jellium_isolated_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/161_Jellium_isolated_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/161_Jellium_isolated_gs/CMakeLists.txt;0;")
add_test(run_161_Jellium_isolated_gs "mpiexec_161_Jellium_isolated_gs.sh")
set_tests_properties(run_161_Jellium_isolated_gs PROPERTIES  FIXTURES_REQUIRED "setup_161_Jellium_isolated_gs" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/161_Jellium_isolated_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/161_Jellium_isolated_gs/CMakeLists.txt;0;")
add_test(verify_161_Jellium_isolated_gs "verification")
set_tests_properties(verify_161_Jellium_isolated_gs PROPERTIES  FIXTURES_CLEANUP "cleanup_161_Jellium_isolated_gs" FIXTURES_REQUIRED "setup_161_Jellium_isolated_gs;run_161_Jellium_isolated_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/161_Jellium_isolated_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/161_Jellium_isolated_gs/CMakeLists.txt;0;")
