# CMake generated Testfile for 
# Source directory: /workspace/testsuites/191_bulk_Si_pseudo_fhi
# Build directory: /workspace/build/testsuites/191_bulk_Si_pseudo_fhi
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_191_bulk_Si_pseudo_fhi "prep_191_bulk_Si_pseudo_fhi.sh")
set_tests_properties(prep_191_bulk_Si_pseudo_fhi PROPERTIES  FIXTURES_SETUP "setup_191_bulk_Si_pseudo_fhi" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/191_bulk_Si_pseudo_fhi/CMakeLists.txt;1;create_test;/workspace/testsuites/191_bulk_Si_pseudo_fhi/CMakeLists.txt;0;")
add_test(run_191_bulk_Si_pseudo_fhi "mpiexec_191_bulk_Si_pseudo_fhi.sh")
set_tests_properties(run_191_bulk_Si_pseudo_fhi PROPERTIES  FIXTURES_REQUIRED "setup_191_bulk_Si_pseudo_fhi" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/191_bulk_Si_pseudo_fhi/CMakeLists.txt;1;create_test;/workspace/testsuites/191_bulk_Si_pseudo_fhi/CMakeLists.txt;0;")
add_test(verify_191_bulk_Si_pseudo_fhi "verification")
set_tests_properties(verify_191_bulk_Si_pseudo_fhi PROPERTIES  FIXTURES_CLEANUP "cleanup_191_bulk_Si_pseudo_fhi" FIXTURES_REQUIRED "setup_191_bulk_Si_pseudo_fhi;run_191_bulk_Si_pseudo_fhi" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/191_bulk_Si_pseudo_fhi/CMakeLists.txt;1;create_test;/workspace/testsuites/191_bulk_Si_pseudo_fhi/CMakeLists.txt;0;")
