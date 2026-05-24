# CMake generated Testfile for 
# Source directory: /workspace/testsuites/103_C2H2_rt_pulse
# Build directory: /workspace/build/testsuites/103_C2H2_rt_pulse
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_103_C2H2_rt_pulse "prep_103_C2H2_rt_pulse.sh")
set_tests_properties(prep_103_C2H2_rt_pulse PROPERTIES  FIXTURES_SETUP "setup_103_C2H2_rt_pulse" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/103_C2H2_rt_pulse/CMakeLists.txt;1;create_test;/workspace/testsuites/103_C2H2_rt_pulse/CMakeLists.txt;0;")
add_test(run_103_C2H2_rt_pulse "mpiexec_103_C2H2_rt_pulse.sh")
set_tests_properties(run_103_C2H2_rt_pulse PROPERTIES  FIXTURES_REQUIRED "setup_103_C2H2_rt_pulse" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/103_C2H2_rt_pulse/CMakeLists.txt;1;create_test;/workspace/testsuites/103_C2H2_rt_pulse/CMakeLists.txt;0;")
add_test(verify_103_C2H2_rt_pulse "verification")
set_tests_properties(verify_103_C2H2_rt_pulse PROPERTIES  FIXTURES_CLEANUP "cleanup_103_C2H2_rt_pulse" FIXTURES_REQUIRED "setup_103_C2H2_rt_pulse;run_103_C2H2_rt_pulse" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/103_C2H2_rt_pulse/CMakeLists.txt;1;create_test;/workspace/testsuites/103_C2H2_rt_pulse/CMakeLists.txt;0;")
