# CMake generated Testfile for 
# Source directory: /workspace/testsuites/141_bulk_Si_singlescale
# Build directory: /workspace/build/testsuites/141_bulk_Si_singlescale
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_141_bulk_Si_singlescale "prep_141_bulk_Si_singlescale.sh")
set_tests_properties(prep_141_bulk_Si_singlescale PROPERTIES  FIXTURES_SETUP "setup_141_bulk_Si_singlescale" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/141_bulk_Si_singlescale/CMakeLists.txt;1;create_test;/workspace/testsuites/141_bulk_Si_singlescale/CMakeLists.txt;0;")
add_test(run_141_bulk_Si_singlescale "mpiexec_141_bulk_Si_singlescale.sh")
set_tests_properties(run_141_bulk_Si_singlescale PROPERTIES  FIXTURES_REQUIRED "setup_141_bulk_Si_singlescale" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/141_bulk_Si_singlescale/CMakeLists.txt;1;create_test;/workspace/testsuites/141_bulk_Si_singlescale/CMakeLists.txt;0;")
add_test(verify_141_bulk_Si_singlescale "verification")
set_tests_properties(verify_141_bulk_Si_singlescale PROPERTIES  FIXTURES_CLEANUP "cleanup_141_bulk_Si_singlescale" FIXTURES_REQUIRED "setup_141_bulk_Si_singlescale;run_141_bulk_Si_singlescale" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/141_bulk_Si_singlescale/CMakeLists.txt;1;create_test;/workspace/testsuites/141_bulk_Si_singlescale/CMakeLists.txt;0;")
