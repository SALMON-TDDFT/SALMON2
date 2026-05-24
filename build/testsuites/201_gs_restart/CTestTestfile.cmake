# CMake generated Testfile for 
# Source directory: /workspace/testsuites/201_gs_restart
# Build directory: /workspace/build/testsuites/201_gs_restart
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_201_gs_restart "prep_201_gs_restart.sh")
set_tests_properties(prep_201_gs_restart PROPERTIES  FIXTURES_SETUP "setup_201_gs_restart" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/201_gs_restart/CMakeLists.txt;1;create_test;/workspace/testsuites/201_gs_restart/CMakeLists.txt;0;")
add_test(run_201_gs_restart "mpiexec_201_gs_restart.sh")
set_tests_properties(run_201_gs_restart PROPERTIES  FIXTURES_REQUIRED "setup_201_gs_restart" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/201_gs_restart/CMakeLists.txt;1;create_test;/workspace/testsuites/201_gs_restart/CMakeLists.txt;0;")
add_test(verify_201_gs_restart "verification")
set_tests_properties(verify_201_gs_restart PROPERTIES  FIXTURES_CLEANUP "cleanup_201_gs_restart" FIXTURES_REQUIRED "setup_201_gs_restart;run_201_gs_restart" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/201_gs_restart/CMakeLists.txt;1;create_test;/workspace/testsuites/201_gs_restart/CMakeLists.txt;0;")
