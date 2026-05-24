# CMake generated Testfile for 
# Source directory: /workspace/testsuites/702_C2H2_rt_response_dirichlet
# Build directory: /workspace/build/testsuites/702_C2H2_rt_response_dirichlet
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_702_C2H2_rt_response_dirichlet "prep_702_C2H2_rt_response_dirichlet.sh")
set_tests_properties(prep_702_C2H2_rt_response_dirichlet PROPERTIES  FIXTURES_SETUP "setup_702_C2H2_rt_response_dirichlet" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/702_C2H2_rt_response_dirichlet/CMakeLists.txt;1;create_test;/workspace/testsuites/702_C2H2_rt_response_dirichlet/CMakeLists.txt;0;")
add_test(run_702_C2H2_rt_response_dirichlet "mpiexec_702_C2H2_rt_response_dirichlet.sh")
set_tests_properties(run_702_C2H2_rt_response_dirichlet PROPERTIES  FIXTURES_REQUIRED "setup_702_C2H2_rt_response_dirichlet" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/702_C2H2_rt_response_dirichlet/CMakeLists.txt;1;create_test;/workspace/testsuites/702_C2H2_rt_response_dirichlet/CMakeLists.txt;0;")
add_test(verify_702_C2H2_rt_response_dirichlet "verification")
set_tests_properties(verify_702_C2H2_rt_response_dirichlet PROPERTIES  FIXTURES_CLEANUP "cleanup_702_C2H2_rt_response_dirichlet" FIXTURES_REQUIRED "setup_702_C2H2_rt_response_dirichlet;run_702_C2H2_rt_response_dirichlet" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/702_C2H2_rt_response_dirichlet/CMakeLists.txt;1;create_test;/workspace/testsuites/702_C2H2_rt_response_dirichlet/CMakeLists.txt;0;")
