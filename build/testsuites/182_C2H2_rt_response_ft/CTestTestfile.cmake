# CMake generated Testfile for 
# Source directory: /workspace/testsuites/182_C2H2_rt_response_ft
# Build directory: /workspace/build/testsuites/182_C2H2_rt_response_ft
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_182_C2H2_rt_response_ft "prep_182_C2H2_rt_response_ft.sh")
set_tests_properties(prep_182_C2H2_rt_response_ft PROPERTIES  FIXTURES_SETUP "setup_182_C2H2_rt_response_ft" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/182_C2H2_rt_response_ft/CMakeLists.txt;1;create_test;/workspace/testsuites/182_C2H2_rt_response_ft/CMakeLists.txt;0;")
add_test(run_182_C2H2_rt_response_ft "mpiexec_182_C2H2_rt_response_ft.sh")
set_tests_properties(run_182_C2H2_rt_response_ft PROPERTIES  FIXTURES_REQUIRED "setup_182_C2H2_rt_response_ft" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/182_C2H2_rt_response_ft/CMakeLists.txt;1;create_test;/workspace/testsuites/182_C2H2_rt_response_ft/CMakeLists.txt;0;")
add_test(verify_182_C2H2_rt_response_ft "verification")
set_tests_properties(verify_182_C2H2_rt_response_ft PROPERTIES  FIXTURES_CLEANUP "cleanup_182_C2H2_rt_response_ft" FIXTURES_REQUIRED "setup_182_C2H2_rt_response_ft;run_182_C2H2_rt_response_ft" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/182_C2H2_rt_response_ft/CMakeLists.txt;1;create_test;/workspace/testsuites/182_C2H2_rt_response_ft/CMakeLists.txt;0;")
