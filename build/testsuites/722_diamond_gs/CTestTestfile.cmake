# CMake generated Testfile for 
# Source directory: /workspace/testsuites/722_diamond_gs
# Build directory: /workspace/build/testsuites/722_diamond_gs
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_722_diamond_gs "prep_722_diamond_gs.sh")
set_tests_properties(prep_722_diamond_gs PROPERTIES  FIXTURES_SETUP "setup_722_diamond_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/722_diamond_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/722_diamond_gs/CMakeLists.txt;0;")
add_test(run_722_diamond_gs "mpiexec_722_diamond_gs.sh")
set_tests_properties(run_722_diamond_gs PROPERTIES  FIXTURES_REQUIRED "setup_722_diamond_gs" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/722_diamond_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/722_diamond_gs/CMakeLists.txt;0;")
add_test(verify_722_diamond_gs "verification")
set_tests_properties(verify_722_diamond_gs PROPERTIES  FIXTURES_CLEANUP "cleanup_722_diamond_gs" FIXTURES_REQUIRED "setup_722_diamond_gs;run_722_diamond_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/722_diamond_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/722_diamond_gs/CMakeLists.txt;0;")
