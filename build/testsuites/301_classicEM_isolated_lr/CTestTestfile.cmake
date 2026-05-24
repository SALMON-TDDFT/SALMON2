# CMake generated Testfile for 
# Source directory: /workspace/testsuites/301_classicEM_isolated_lr
# Build directory: /workspace/build/testsuites/301_classicEM_isolated_lr
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_301_classicEM_isolated_lr "prep_301_classicEM_isolated_lr.sh")
set_tests_properties(prep_301_classicEM_isolated_lr PROPERTIES  FIXTURES_SETUP "setup_301_classicEM_isolated_lr" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/301_classicEM_isolated_lr/CMakeLists.txt;1;create_test;/workspace/testsuites/301_classicEM_isolated_lr/CMakeLists.txt;0;")
add_test(run_301_classicEM_isolated_lr "mpiexec_301_classicEM_isolated_lr.sh")
set_tests_properties(run_301_classicEM_isolated_lr PROPERTIES  FIXTURES_REQUIRED "setup_301_classicEM_isolated_lr" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/301_classicEM_isolated_lr/CMakeLists.txt;1;create_test;/workspace/testsuites/301_classicEM_isolated_lr/CMakeLists.txt;0;")
add_test(verify_301_classicEM_isolated_lr "verification")
set_tests_properties(verify_301_classicEM_isolated_lr PROPERTIES  FIXTURES_CLEANUP "cleanup_301_classicEM_isolated_lr" FIXTURES_REQUIRED "setup_301_classicEM_isolated_lr;run_301_classicEM_isolated_lr" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/301_classicEM_isolated_lr/CMakeLists.txt;1;create_test;/workspace/testsuites/301_classicEM_isolated_lr/CMakeLists.txt;0;")
