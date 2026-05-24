# CMake generated Testfile for 
# Source directory: /workspace/testsuites/410_bulk_Si_gs_initial_density
# Build directory: /workspace/build/testsuites/410_bulk_Si_gs_initial_density
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_410_bulk_Si_gs_initial_density "prep_410_bulk_Si_gs_initial_density.sh")
set_tests_properties(prep_410_bulk_Si_gs_initial_density PROPERTIES  FIXTURES_SETUP "setup_410_bulk_Si_gs_initial_density" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/410_bulk_Si_gs_initial_density/CMakeLists.txt;1;create_test;/workspace/testsuites/410_bulk_Si_gs_initial_density/CMakeLists.txt;0;")
add_test(run_410_bulk_Si_gs_initial_density "mpiexec_410_bulk_Si_gs_initial_density.sh")
set_tests_properties(run_410_bulk_Si_gs_initial_density PROPERTIES  FIXTURES_REQUIRED "setup_410_bulk_Si_gs_initial_density" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/410_bulk_Si_gs_initial_density/CMakeLists.txt;1;create_test;/workspace/testsuites/410_bulk_Si_gs_initial_density/CMakeLists.txt;0;")
add_test(verify_410_bulk_Si_gs_initial_density "verification")
set_tests_properties(verify_410_bulk_Si_gs_initial_density PROPERTIES  FIXTURES_CLEANUP "cleanup_410_bulk_Si_gs_initial_density" FIXTURES_REQUIRED "setup_410_bulk_Si_gs_initial_density;run_410_bulk_Si_gs_initial_density" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/410_bulk_Si_gs_initial_density/CMakeLists.txt;1;create_test;/workspace/testsuites/410_bulk_Si_gs_initial_density/CMakeLists.txt;0;")
