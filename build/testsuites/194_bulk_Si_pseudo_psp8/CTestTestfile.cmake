# CMake generated Testfile for 
# Source directory: /workspace/testsuites/194_bulk_Si_pseudo_psp8
# Build directory: /workspace/build/testsuites/194_bulk_Si_pseudo_psp8
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_194_bulk_Si_pseudo_psp8 "prep_194_bulk_Si_pseudo_psp8.sh")
set_tests_properties(prep_194_bulk_Si_pseudo_psp8 PROPERTIES  FIXTURES_SETUP "setup_194_bulk_Si_pseudo_psp8" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/194_bulk_Si_pseudo_psp8/CMakeLists.txt;1;create_test;/workspace/testsuites/194_bulk_Si_pseudo_psp8/CMakeLists.txt;0;")
add_test(run_194_bulk_Si_pseudo_psp8 "mpiexec_194_bulk_Si_pseudo_psp8.sh")
set_tests_properties(run_194_bulk_Si_pseudo_psp8 PROPERTIES  FIXTURES_REQUIRED "setup_194_bulk_Si_pseudo_psp8" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/194_bulk_Si_pseudo_psp8/CMakeLists.txt;1;create_test;/workspace/testsuites/194_bulk_Si_pseudo_psp8/CMakeLists.txt;0;")
add_test(verify_194_bulk_Si_pseudo_psp8 "verification")
set_tests_properties(verify_194_bulk_Si_pseudo_psp8 PROPERTIES  FIXTURES_CLEANUP "cleanup_194_bulk_Si_pseudo_psp8" FIXTURES_REQUIRED "setup_194_bulk_Si_pseudo_psp8;run_194_bulk_Si_pseudo_psp8" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/194_bulk_Si_pseudo_psp8/CMakeLists.txt;1;create_test;/workspace/testsuites/194_bulk_Si_pseudo_psp8/CMakeLists.txt;0;")
