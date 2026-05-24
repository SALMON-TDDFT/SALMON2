# CMake generated Testfile for 
# Source directory: /workspace/testsuites/711_bulk_gs_pre_dp
# Build directory: /workspace/build/testsuites/711_bulk_gs_pre_dp
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_711_bulk_gs_pre_dp "prep_711_bulk_gs_pre_dp.sh")
set_tests_properties(prep_711_bulk_gs_pre_dp PROPERTIES  FIXTURES_SETUP "setup_711_bulk_gs_pre_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/711_bulk_gs_pre_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/711_bulk_gs_pre_dp/CMakeLists.txt;0;")
add_test(run_711_bulk_gs_pre_dp "mpiexec_711_bulk_gs_pre_dp.sh")
set_tests_properties(run_711_bulk_gs_pre_dp PROPERTIES  FIXTURES_REQUIRED "setup_711_bulk_gs_pre_dp" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/711_bulk_gs_pre_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/711_bulk_gs_pre_dp/CMakeLists.txt;0;")
add_test(verify_711_bulk_gs_pre_dp "verification")
set_tests_properties(verify_711_bulk_gs_pre_dp PROPERTIES  FIXTURES_CLEANUP "cleanup_711_bulk_gs_pre_dp" FIXTURES_REQUIRED "setup_711_bulk_gs_pre_dp;run_711_bulk_gs_pre_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/711_bulk_gs_pre_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/711_bulk_gs_pre_dp/CMakeLists.txt;0;")
