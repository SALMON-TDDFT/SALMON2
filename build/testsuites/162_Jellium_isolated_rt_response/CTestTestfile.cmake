# CMake generated Testfile for 
# Source directory: /workspace/testsuites/162_Jellium_isolated_rt_response
# Build directory: /workspace/build/testsuites/162_Jellium_isolated_rt_response
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_162_Jellium_isolated_rt_response "prep_162_Jellium_isolated_rt_response.sh")
set_tests_properties(prep_162_Jellium_isolated_rt_response PROPERTIES  FIXTURES_SETUP "setup_162_Jellium_isolated_rt_response" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/162_Jellium_isolated_rt_response/CMakeLists.txt;1;create_test;/workspace/testsuites/162_Jellium_isolated_rt_response/CMakeLists.txt;0;")
add_test(run_162_Jellium_isolated_rt_response "mpiexec_162_Jellium_isolated_rt_response.sh")
set_tests_properties(run_162_Jellium_isolated_rt_response PROPERTIES  FIXTURES_REQUIRED "setup_162_Jellium_isolated_rt_response" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/162_Jellium_isolated_rt_response/CMakeLists.txt;1;create_test;/workspace/testsuites/162_Jellium_isolated_rt_response/CMakeLists.txt;0;")
add_test(verify_162_Jellium_isolated_rt_response "verification")
set_tests_properties(verify_162_Jellium_isolated_rt_response PROPERTIES  FIXTURES_CLEANUP "cleanup_162_Jellium_isolated_rt_response" FIXTURES_REQUIRED "setup_162_Jellium_isolated_rt_response;run_162_Jellium_isolated_rt_response" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/162_Jellium_isolated_rt_response/CMakeLists.txt;1;create_test;/workspace/testsuites/162_Jellium_isolated_rt_response/CMakeLists.txt;0;")
