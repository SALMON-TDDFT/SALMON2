# CMake generated Testfile for 
# Source directory: /workspace/testsuites/202_rt_checkpoint_single
# Build directory: /workspace/build/testsuites/202_rt_checkpoint_single
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_202_rt_checkpoint_single "prep_202_rt_checkpoint_single.sh")
set_tests_properties(prep_202_rt_checkpoint_single PROPERTIES  FIXTURES_SETUP "setup_202_rt_checkpoint_single" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/202_rt_checkpoint_single/CMakeLists.txt;1;create_test;/workspace/testsuites/202_rt_checkpoint_single/CMakeLists.txt;0;")
add_test(run_202_rt_checkpoint_single "mpiexec_202_rt_checkpoint_single.sh")
set_tests_properties(run_202_rt_checkpoint_single PROPERTIES  FIXTURES_REQUIRED "setup_202_rt_checkpoint_single" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/202_rt_checkpoint_single/CMakeLists.txt;1;create_test;/workspace/testsuites/202_rt_checkpoint_single/CMakeLists.txt;0;")
add_test(verify_202_rt_checkpoint_single "verification")
set_tests_properties(verify_202_rt_checkpoint_single PROPERTIES  FIXTURES_CLEANUP "cleanup_202_rt_checkpoint_single" FIXTURES_REQUIRED "setup_202_rt_checkpoint_single;run_202_rt_checkpoint_single" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/202_rt_checkpoint_single/CMakeLists.txt;1;create_test;/workspace/testsuites/202_rt_checkpoint_single/CMakeLists.txt;0;")
