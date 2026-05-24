# CMake generated Testfile for 
# Source directory: /workspace/testsuites/502_bulk_Au_spin-orbit_rt
# Build directory: /workspace/build/testsuites/502_bulk_Au_spin-orbit_rt
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_502_bulk_Au_spin-orbit_rt "prep_502_bulk_Au_spin-orbit_rt.sh")
set_tests_properties(prep_502_bulk_Au_spin-orbit_rt PROPERTIES  FIXTURES_SETUP "setup_502_bulk_Au_spin-orbit_rt" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/502_bulk_Au_spin-orbit_rt/CMakeLists.txt;1;create_test;/workspace/testsuites/502_bulk_Au_spin-orbit_rt/CMakeLists.txt;0;")
add_test(run_502_bulk_Au_spin-orbit_rt "mpiexec_502_bulk_Au_spin-orbit_rt.sh")
set_tests_properties(run_502_bulk_Au_spin-orbit_rt PROPERTIES  FIXTURES_REQUIRED "setup_502_bulk_Au_spin-orbit_rt" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/502_bulk_Au_spin-orbit_rt/CMakeLists.txt;1;create_test;/workspace/testsuites/502_bulk_Au_spin-orbit_rt/CMakeLists.txt;0;")
add_test(verify_502_bulk_Au_spin-orbit_rt "verification")
set_tests_properties(verify_502_bulk_Au_spin-orbit_rt PROPERTIES  FIXTURES_CLEANUP "cleanup_502_bulk_Au_spin-orbit_rt" FIXTURES_REQUIRED "setup_502_bulk_Au_spin-orbit_rt;run_502_bulk_Au_spin-orbit_rt" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/502_bulk_Au_spin-orbit_rt/CMakeLists.txt;1;create_test;/workspace/testsuites/502_bulk_Au_spin-orbit_rt/CMakeLists.txt;0;")
