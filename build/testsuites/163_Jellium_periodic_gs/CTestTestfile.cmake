# CMake generated Testfile for 
# Source directory: /workspace/testsuites/163_Jellium_periodic_gs
# Build directory: /workspace/build/testsuites/163_Jellium_periodic_gs
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_163_Jellium_periodic_gs "prep_163_Jellium_periodic_gs.sh")
set_tests_properties(prep_163_Jellium_periodic_gs PROPERTIES  FIXTURES_SETUP "setup_163_Jellium_periodic_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/163_Jellium_periodic_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/163_Jellium_periodic_gs/CMakeLists.txt;0;")
add_test(run_163_Jellium_periodic_gs "mpiexec_163_Jellium_periodic_gs.sh")
set_tests_properties(run_163_Jellium_periodic_gs PROPERTIES  FIXTURES_REQUIRED "setup_163_Jellium_periodic_gs" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/163_Jellium_periodic_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/163_Jellium_periodic_gs/CMakeLists.txt;0;")
add_test(verify_163_Jellium_periodic_gs "verification")
set_tests_properties(verify_163_Jellium_periodic_gs PROPERTIES  FIXTURES_CLEANUP "cleanup_163_Jellium_periodic_gs" FIXTURES_REQUIRED "setup_163_Jellium_periodic_gs;run_163_Jellium_periodic_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/163_Jellium_periodic_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/163_Jellium_periodic_gs/CMakeLists.txt;0;")
