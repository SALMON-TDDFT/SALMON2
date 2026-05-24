# CMake generated Testfile for 
# Source directory: /workspace/testsuites/101_C2H2_gs
# Build directory: /workspace/build/testsuites/101_C2H2_gs
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_101_C2H2_gs "prep_101_C2H2_gs.sh")
set_tests_properties(prep_101_C2H2_gs PROPERTIES  FIXTURES_SETUP "setup_101_C2H2_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/101_C2H2_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/101_C2H2_gs/CMakeLists.txt;0;")
add_test(run_101_C2H2_gs "mpiexec_101_C2H2_gs.sh")
set_tests_properties(run_101_C2H2_gs PROPERTIES  FIXTURES_REQUIRED "setup_101_C2H2_gs" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/101_C2H2_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/101_C2H2_gs/CMakeLists.txt;0;")
add_test(verify_101_C2H2_gs "verification")
set_tests_properties(verify_101_C2H2_gs PROPERTIES  FIXTURES_CLEANUP "cleanup_101_C2H2_gs" FIXTURES_REQUIRED "setup_101_C2H2_gs;run_101_C2H2_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/101_C2H2_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/101_C2H2_gs/CMakeLists.txt;0;")
