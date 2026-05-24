# CMake generated Testfile for 
# Source directory: /workspace/testsuites/207_rt_wf_slice
# Build directory: /workspace/build/testsuites/207_rt_wf_slice
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_207_rt_wf_slice "prep_207_rt_wf_slice.sh")
set_tests_properties(prep_207_rt_wf_slice PROPERTIES  FIXTURES_SETUP "setup_207_rt_wf_slice" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/207_rt_wf_slice/CMakeLists.txt;1;create_test;/workspace/testsuites/207_rt_wf_slice/CMakeLists.txt;0;")
add_test(run_207_rt_wf_slice "mpiexec_207_rt_wf_slice.sh")
set_tests_properties(run_207_rt_wf_slice PROPERTIES  FIXTURES_REQUIRED "setup_207_rt_wf_slice" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/207_rt_wf_slice/CMakeLists.txt;1;create_test;/workspace/testsuites/207_rt_wf_slice/CMakeLists.txt;0;")
add_test(verify_207_rt_wf_slice "verification")
set_tests_properties(verify_207_rt_wf_slice PROPERTIES  FIXTURES_CLEANUP "cleanup_207_rt_wf_slice" FIXTURES_REQUIRED "setup_207_rt_wf_slice;run_207_rt_wf_slice" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/207_rt_wf_slice/CMakeLists.txt;1;create_test;/workspace/testsuites/207_rt_wf_slice/CMakeLists.txt;0;")
