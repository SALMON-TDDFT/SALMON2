# CMake generated Testfile for 
# Source directory: /workspace/testsuites/203_rt_restart_single
# Build directory: /workspace/build/testsuites/203_rt_restart_single
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_203_rt_restart_single "prep_203_rt_restart_single.sh")
set_tests_properties(prep_203_rt_restart_single PROPERTIES  FIXTURES_SETUP "setup_203_rt_restart_single" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/203_rt_restart_single/CMakeLists.txt;1;create_test;/workspace/testsuites/203_rt_restart_single/CMakeLists.txt;0;")
add_test(run_203_rt_restart_single "mpiexec_203_rt_restart_single.sh")
set_tests_properties(run_203_rt_restart_single PROPERTIES  FIXTURES_REQUIRED "setup_203_rt_restart_single" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/203_rt_restart_single/CMakeLists.txt;1;create_test;/workspace/testsuites/203_rt_restart_single/CMakeLists.txt;0;")
add_test(verify_203_rt_restart_single "verification")
set_tests_properties(verify_203_rt_restart_single PROPERTIES  FIXTURES_CLEANUP "cleanup_203_rt_restart_single" FIXTURES_REQUIRED "setup_203_rt_restart_single;run_203_rt_restart_single" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/203_rt_restart_single/CMakeLists.txt;1;create_test;/workspace/testsuites/203_rt_restart_single/CMakeLists.txt;0;")
