# CMake generated Testfile for 
# Source directory: /workspace/testsuites/124_periodic_H2O_md_dp
# Build directory: /workspace/build/testsuites/124_periodic_H2O_md_dp
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_124_periodic_H2O_md_dp "prep_124_periodic_H2O_md_dp.sh")
set_tests_properties(prep_124_periodic_H2O_md_dp PROPERTIES  FIXTURES_SETUP "setup_124_periodic_H2O_md_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/124_periodic_H2O_md_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/124_periodic_H2O_md_dp/CMakeLists.txt;0;")
add_test(run_124_periodic_H2O_md_dp "mpiexec_124_periodic_H2O_md_dp.sh")
set_tests_properties(run_124_periodic_H2O_md_dp PROPERTIES  FIXTURES_REQUIRED "setup_124_periodic_H2O_md_dp" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/124_periodic_H2O_md_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/124_periodic_H2O_md_dp/CMakeLists.txt;0;")
add_test(verify_124_periodic_H2O_md_dp "verification")
set_tests_properties(verify_124_periodic_H2O_md_dp PROPERTIES  FIXTURES_CLEANUP "cleanup_124_periodic_H2O_md_dp" FIXTURES_REQUIRED "setup_124_periodic_H2O_md_dp;run_124_periodic_H2O_md_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/124_periodic_H2O_md_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/124_periodic_H2O_md_dp/CMakeLists.txt;0;")
