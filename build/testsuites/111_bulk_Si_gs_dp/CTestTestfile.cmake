# CMake generated Testfile for 
# Source directory: /workspace/testsuites/111_bulk_Si_gs_dp
# Build directory: /workspace/build/testsuites/111_bulk_Si_gs_dp
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_111_bulk_Si_gs_dp "prep_111_bulk_Si_gs_dp.sh")
set_tests_properties(prep_111_bulk_Si_gs_dp PROPERTIES  FIXTURES_SETUP "setup_111_bulk_Si_gs_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/111_bulk_Si_gs_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/111_bulk_Si_gs_dp/CMakeLists.txt;0;")
add_test(run_111_bulk_Si_gs_dp "mpiexec_111_bulk_Si_gs_dp.sh")
set_tests_properties(run_111_bulk_Si_gs_dp PROPERTIES  FIXTURES_REQUIRED "setup_111_bulk_Si_gs_dp" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/111_bulk_Si_gs_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/111_bulk_Si_gs_dp/CMakeLists.txt;0;")
add_test(verify_111_bulk_Si_gs_dp "verification")
set_tests_properties(verify_111_bulk_Si_gs_dp PROPERTIES  FIXTURES_CLEANUP "cleanup_111_bulk_Si_gs_dp" FIXTURES_REQUIRED "setup_111_bulk_Si_gs_dp;run_111_bulk_Si_gs_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/111_bulk_Si_gs_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/111_bulk_Si_gs_dp/CMakeLists.txt;0;")
