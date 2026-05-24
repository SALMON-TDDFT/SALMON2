# CMake generated Testfile for 
# Source directory: /workspace/testsuites/119_bulk_Si_rt_response_dp_ffte
# Build directory: /workspace/build/testsuites/119_bulk_Si_rt_response_dp_ffte
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_119_bulk_Si_rt_response_dp_ffte "prep_119_bulk_Si_rt_response_dp_ffte.sh")
set_tests_properties(prep_119_bulk_Si_rt_response_dp_ffte PROPERTIES  FIXTURES_SETUP "setup_119_bulk_Si_rt_response_dp_ffte" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/119_bulk_Si_rt_response_dp_ffte/CMakeLists.txt;1;create_test;/workspace/testsuites/119_bulk_Si_rt_response_dp_ffte/CMakeLists.txt;0;")
add_test(run_119_bulk_Si_rt_response_dp_ffte "mpiexec_119_bulk_Si_rt_response_dp_ffte.sh")
set_tests_properties(run_119_bulk_Si_rt_response_dp_ffte PROPERTIES  FIXTURES_REQUIRED "setup_119_bulk_Si_rt_response_dp_ffte" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/119_bulk_Si_rt_response_dp_ffte/CMakeLists.txt;1;create_test;/workspace/testsuites/119_bulk_Si_rt_response_dp_ffte/CMakeLists.txt;0;")
add_test(verify_119_bulk_Si_rt_response_dp_ffte "verification")
set_tests_properties(verify_119_bulk_Si_rt_response_dp_ffte PROPERTIES  FIXTURES_CLEANUP "cleanup_119_bulk_Si_rt_response_dp_ffte" FIXTURES_REQUIRED "setup_119_bulk_Si_rt_response_dp_ffte;run_119_bulk_Si_rt_response_dp_ffte" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/119_bulk_Si_rt_response_dp_ffte/CMakeLists.txt;1;create_test;/workspace/testsuites/119_bulk_Si_rt_response_dp_ffte/CMakeLists.txt;0;")
