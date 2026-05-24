# CMake generated Testfile for 
# Source directory: /workspace/testsuites/113_bulk_Si_rt_pulse_dp
# Build directory: /workspace/build/testsuites/113_bulk_Si_rt_pulse_dp
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_113_bulk_Si_rt_pulse_dp "prep_113_bulk_Si_rt_pulse_dp.sh")
set_tests_properties(prep_113_bulk_Si_rt_pulse_dp PROPERTIES  FIXTURES_SETUP "setup_113_bulk_Si_rt_pulse_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/113_bulk_Si_rt_pulse_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/113_bulk_Si_rt_pulse_dp/CMakeLists.txt;0;")
add_test(run_113_bulk_Si_rt_pulse_dp "mpiexec_113_bulk_Si_rt_pulse_dp.sh")
set_tests_properties(run_113_bulk_Si_rt_pulse_dp PROPERTIES  FIXTURES_REQUIRED "setup_113_bulk_Si_rt_pulse_dp" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/113_bulk_Si_rt_pulse_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/113_bulk_Si_rt_pulse_dp/CMakeLists.txt;0;")
add_test(verify_113_bulk_Si_rt_pulse_dp "verification")
set_tests_properties(verify_113_bulk_Si_rt_pulse_dp PROPERTIES  FIXTURES_CLEANUP "cleanup_113_bulk_Si_rt_pulse_dp" FIXTURES_REQUIRED "setup_113_bulk_Si_rt_pulse_dp;run_113_bulk_Si_rt_pulse_dp" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/113_bulk_Si_rt_pulse_dp/CMakeLists.txt;1;create_test;/workspace/testsuites/113_bulk_Si_rt_pulse_dp/CMakeLists.txt;0;")
