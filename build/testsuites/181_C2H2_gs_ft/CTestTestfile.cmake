# CMake generated Testfile for 
# Source directory: /workspace/testsuites/181_C2H2_gs_ft
# Build directory: /workspace/build/testsuites/181_C2H2_gs_ft
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_181_C2H2_gs_ft "prep_181_C2H2_gs_ft.sh")
set_tests_properties(prep_181_C2H2_gs_ft PROPERTIES  FIXTURES_SETUP "setup_181_C2H2_gs_ft" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/181_C2H2_gs_ft/CMakeLists.txt;1;create_test;/workspace/testsuites/181_C2H2_gs_ft/CMakeLists.txt;0;")
add_test(run_181_C2H2_gs_ft "mpiexec_181_C2H2_gs_ft.sh")
set_tests_properties(run_181_C2H2_gs_ft PROPERTIES  FIXTURES_REQUIRED "setup_181_C2H2_gs_ft" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/181_C2H2_gs_ft/CMakeLists.txt;1;create_test;/workspace/testsuites/181_C2H2_gs_ft/CMakeLists.txt;0;")
add_test(verify_181_C2H2_gs_ft "verification")
set_tests_properties(verify_181_C2H2_gs_ft PROPERTIES  FIXTURES_CLEANUP "cleanup_181_C2H2_gs_ft" FIXTURES_REQUIRED "setup_181_C2H2_gs_ft;run_181_C2H2_gs_ft" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/181_C2H2_gs_ft/CMakeLists.txt;1;create_test;/workspace/testsuites/181_C2H2_gs_ft/CMakeLists.txt;0;")
