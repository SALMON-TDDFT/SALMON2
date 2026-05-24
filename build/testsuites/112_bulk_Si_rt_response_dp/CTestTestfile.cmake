# CMake generated Testfile for 
# Source directory: /workspace/testsuites/112_bulk_Si_rt_response_dp
# Build directory: /workspace/build/testsuites/112_bulk_Si_rt_response_dp
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_112_bulk_Si_rt_response_dp "prep_112_bulk_Si_rt_response_dp.sh")
set_tests_properties(prep_112_bulk_Si_rt_response_dp PROPERTIES  FIXTURES_SETUP "setup_112_bulk_Si_rt_response_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/112_bulk_Si_rt_response_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/112_bulk_Si_rt_response_dp/CMakeLists.txt;0;")
add_test(run_112_bulk_Si_rt_response_dp "mpiexec_112_bulk_Si_rt_response_dp.sh")
set_tests_properties(run_112_bulk_Si_rt_response_dp PROPERTIES  FIXTURES_REQUIRED "setup_112_bulk_Si_rt_response_dp" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/112_bulk_Si_rt_response_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/112_bulk_Si_rt_response_dp/CMakeLists.txt;0;")
add_test(verify_112_bulk_Si_rt_response_dp "verification")
set_tests_properties(verify_112_bulk_Si_rt_response_dp PROPERTIES  FIXTURES_CLEANUP "cleanup_112_bulk_Si_rt_response_dp" FIXTURES_REQUIRED "setup_112_bulk_Si_rt_response_dp;run_112_bulk_Si_rt_response_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/112_bulk_Si_rt_response_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/112_bulk_Si_rt_response_dp/CMakeLists.txt;0;")
