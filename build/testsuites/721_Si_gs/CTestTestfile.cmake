# CMake generated Testfile for 
# Source directory: /workspace/testsuites/721_Si_gs
# Build directory: /workspace/build/testsuites/721_Si_gs
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_721_Si_gs "prep_721_Si_gs.sh")
set_tests_properties(prep_721_Si_gs PROPERTIES  FIXTURES_SETUP "setup_721_Si_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/721_Si_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/721_Si_gs/CMakeLists.txt;0;")
add_test(run_721_Si_gs "mpiexec_721_Si_gs.sh")
set_tests_properties(run_721_Si_gs PROPERTIES  FIXTURES_REQUIRED "setup_721_Si_gs" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/721_Si_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/721_Si_gs/CMakeLists.txt;0;")
add_test(verify_721_Si_gs "verification")
set_tests_properties(verify_721_Si_gs PROPERTIES  FIXTURES_CLEANUP "cleanup_721_Si_gs" FIXTURES_REQUIRED "setup_721_Si_gs;run_721_Si_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/721_Si_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/721_Si_gs/CMakeLists.txt;0;")
