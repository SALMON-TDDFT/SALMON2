# CMake generated Testfile for 
# Source directory: /workspace/testsuites/511_isolated_Au_spin-orbit_gs
# Build directory: /workspace/build/testsuites/511_isolated_Au_spin-orbit_gs
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_511_isolated_Au_spin-orbit_gs "prep_511_isolated_Au_spin-orbit_gs.sh")
set_tests_properties(prep_511_isolated_Au_spin-orbit_gs PROPERTIES  FIXTURES_SETUP "setup_511_isolated_Au_spin-orbit_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/511_isolated_Au_spin-orbit_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/511_isolated_Au_spin-orbit_gs/CMakeLists.txt;0;")
add_test(run_511_isolated_Au_spin-orbit_gs "mpiexec_511_isolated_Au_spin-orbit_gs.sh")
set_tests_properties(run_511_isolated_Au_spin-orbit_gs PROPERTIES  FIXTURES_REQUIRED "setup_511_isolated_Au_spin-orbit_gs" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/511_isolated_Au_spin-orbit_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/511_isolated_Au_spin-orbit_gs/CMakeLists.txt;0;")
add_test(verify_511_isolated_Au_spin-orbit_gs "verification")
set_tests_properties(verify_511_isolated_Au_spin-orbit_gs PROPERTIES  FIXTURES_CLEANUP "cleanup_511_isolated_Au_spin-orbit_gs" FIXTURES_REQUIRED "setup_511_isolated_Au_spin-orbit_gs;run_511_isolated_Au_spin-orbit_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/511_isolated_Au_spin-orbit_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/511_isolated_Au_spin-orbit_gs/CMakeLists.txt;0;")
