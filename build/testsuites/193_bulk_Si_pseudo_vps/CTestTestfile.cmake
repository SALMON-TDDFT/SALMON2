# CMake generated Testfile for 
# Source directory: /workspace/testsuites/193_bulk_Si_pseudo_vps
# Build directory: /workspace/build/testsuites/193_bulk_Si_pseudo_vps
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_193_bulk_Si_pseudo_vps "prep_193_bulk_Si_pseudo_vps.sh")
set_tests_properties(prep_193_bulk_Si_pseudo_vps PROPERTIES  FIXTURES_SETUP "setup_193_bulk_Si_pseudo_vps" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/193_bulk_Si_pseudo_vps/CMakeLists.txt;1;create_test;/workspace/testsuites/193_bulk_Si_pseudo_vps/CMakeLists.txt;0;")
add_test(run_193_bulk_Si_pseudo_vps "mpiexec_193_bulk_Si_pseudo_vps.sh")
set_tests_properties(run_193_bulk_Si_pseudo_vps PROPERTIES  FIXTURES_REQUIRED "setup_193_bulk_Si_pseudo_vps" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/193_bulk_Si_pseudo_vps/CMakeLists.txt;1;create_test;/workspace/testsuites/193_bulk_Si_pseudo_vps/CMakeLists.txt;0;")
add_test(verify_193_bulk_Si_pseudo_vps "verification")
set_tests_properties(verify_193_bulk_Si_pseudo_vps PROPERTIES  FIXTURES_CLEANUP "cleanup_193_bulk_Si_pseudo_vps" FIXTURES_REQUIRED "setup_193_bulk_Si_pseudo_vps;run_193_bulk_Si_pseudo_vps" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/193_bulk_Si_pseudo_vps/CMakeLists.txt;1;create_test;/workspace/testsuites/193_bulk_Si_pseudo_vps/CMakeLists.txt;0;")
