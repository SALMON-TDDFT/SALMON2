# CMake generated Testfile for 
# Source directory: /workspace/testsuites/501_bulk_Au_spin-orbit_gs
# Build directory: /workspace/build/testsuites/501_bulk_Au_spin-orbit_gs
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_501_bulk_Au_spin-orbit_gs "prep_501_bulk_Au_spin-orbit_gs.sh")
set_tests_properties(prep_501_bulk_Au_spin-orbit_gs PROPERTIES  FIXTURES_SETUP "setup_501_bulk_Au_spin-orbit_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/501_bulk_Au_spin-orbit_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/501_bulk_Au_spin-orbit_gs/CMakeLists.txt;0;")
add_test(run_501_bulk_Au_spin-orbit_gs "mpiexec_501_bulk_Au_spin-orbit_gs.sh")
set_tests_properties(run_501_bulk_Au_spin-orbit_gs PROPERTIES  FIXTURES_REQUIRED "setup_501_bulk_Au_spin-orbit_gs" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/501_bulk_Au_spin-orbit_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/501_bulk_Au_spin-orbit_gs/CMakeLists.txt;0;")
add_test(verify_501_bulk_Au_spin-orbit_gs "verification")
set_tests_properties(verify_501_bulk_Au_spin-orbit_gs PROPERTIES  FIXTURES_CLEANUP "cleanup_501_bulk_Au_spin-orbit_gs" FIXTURES_REQUIRED "setup_501_bulk_Au_spin-orbit_gs;run_501_bulk_Au_spin-orbit_gs" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/501_bulk_Au_spin-orbit_gs/CMakeLists.txt;1;create_test;/workspace/testsuites/501_bulk_Au_spin-orbit_gs/CMakeLists.txt;0;")
