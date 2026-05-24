# CMake generated Testfile for 
# Source directory: /workspace/testsuites/701_C2H2_gs_dirichlet
# Build directory: /workspace/build/testsuites/701_C2H2_gs_dirichlet
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_701_C2H2_gs_dirichlet "prep_701_C2H2_gs_dirichlet.sh")
set_tests_properties(prep_701_C2H2_gs_dirichlet PROPERTIES  FIXTURES_SETUP "setup_701_C2H2_gs_dirichlet" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/701_C2H2_gs_dirichlet/CMakeLists.txt;1;create_test;/workspace/testsuites/701_C2H2_gs_dirichlet/CMakeLists.txt;0;")
add_test(run_701_C2H2_gs_dirichlet "mpiexec_701_C2H2_gs_dirichlet.sh")
set_tests_properties(run_701_C2H2_gs_dirichlet PROPERTIES  FIXTURES_REQUIRED "setup_701_C2H2_gs_dirichlet" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/701_C2H2_gs_dirichlet/CMakeLists.txt;1;create_test;/workspace/testsuites/701_C2H2_gs_dirichlet/CMakeLists.txt;0;")
add_test(verify_701_C2H2_gs_dirichlet "verification")
set_tests_properties(verify_701_C2H2_gs_dirichlet PROPERTIES  FIXTURES_CLEANUP "cleanup_701_C2H2_gs_dirichlet" FIXTURES_REQUIRED "setup_701_C2H2_gs_dirichlet;run_701_C2H2_gs_dirichlet" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/701_C2H2_gs_dirichlet/CMakeLists.txt;1;create_test;/workspace/testsuites/701_C2H2_gs_dirichlet/CMakeLists.txt;0;")
