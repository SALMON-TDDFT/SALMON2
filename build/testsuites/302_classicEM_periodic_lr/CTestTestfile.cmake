# CMake generated Testfile for 
# Source directory: /workspace/testsuites/302_classicEM_periodic_lr
# Build directory: /workspace/build/testsuites/302_classicEM_periodic_lr
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_302_classicEM_periodic_lr "prep_302_classicEM_periodic_lr.sh")
set_tests_properties(prep_302_classicEM_periodic_lr PROPERTIES  FIXTURES_SETUP "setup_302_classicEM_periodic_lr" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/302_classicEM_periodic_lr/CMakeLists.txt;1;create_test;/workspace/testsuites/302_classicEM_periodic_lr/CMakeLists.txt;0;")
add_test(run_302_classicEM_periodic_lr "mpiexec_302_classicEM_periodic_lr.sh")
set_tests_properties(run_302_classicEM_periodic_lr PROPERTIES  FIXTURES_REQUIRED "setup_302_classicEM_periodic_lr" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/302_classicEM_periodic_lr/CMakeLists.txt;1;create_test;/workspace/testsuites/302_classicEM_periodic_lr/CMakeLists.txt;0;")
add_test(verify_302_classicEM_periodic_lr "verification")
set_tests_properties(verify_302_classicEM_periodic_lr PROPERTIES  FIXTURES_CLEANUP "cleanup_302_classicEM_periodic_lr" FIXTURES_REQUIRED "setup_302_classicEM_periodic_lr;run_302_classicEM_periodic_lr" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/302_classicEM_periodic_lr/CMakeLists.txt;1;create_test;/workspace/testsuites/302_classicEM_periodic_lr/CMakeLists.txt;0;")
