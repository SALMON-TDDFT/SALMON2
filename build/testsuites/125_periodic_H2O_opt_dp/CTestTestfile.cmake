# CMake generated Testfile for 
# Source directory: /workspace/testsuites/125_periodic_H2O_opt_dp
# Build directory: /workspace/build/testsuites/125_periodic_H2O_opt_dp
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_125_periodic_H2O_opt_dp "prep_125_periodic_H2O_opt_dp.sh")
set_tests_properties(prep_125_periodic_H2O_opt_dp PROPERTIES  FIXTURES_SETUP "setup_125_periodic_H2O_opt_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/125_periodic_H2O_opt_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/125_periodic_H2O_opt_dp/CMakeLists.txt;0;")
add_test(run_125_periodic_H2O_opt_dp "mpiexec_125_periodic_H2O_opt_dp.sh")
set_tests_properties(run_125_periodic_H2O_opt_dp PROPERTIES  FIXTURES_REQUIRED "setup_125_periodic_H2O_opt_dp" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/125_periodic_H2O_opt_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/125_periodic_H2O_opt_dp/CMakeLists.txt;0;")
add_test(verify_125_periodic_H2O_opt_dp "verification")
set_tests_properties(verify_125_periodic_H2O_opt_dp PROPERTIES  FIXTURES_CLEANUP "cleanup_125_periodic_H2O_opt_dp" FIXTURES_REQUIRED "setup_125_periodic_H2O_opt_dp;run_125_periodic_H2O_opt_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/125_periodic_H2O_opt_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/125_periodic_H2O_opt_dp/CMakeLists.txt;0;")
