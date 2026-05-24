# CMake generated Testfile for 
# Source directory: /workspace/testsuites/411_bulk_Si_gs_symmetry
# Build directory: /workspace/build/testsuites/411_bulk_Si_gs_symmetry
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_411_bulk_Si_gs_symmetry "prep_411_bulk_Si_gs_symmetry.sh")
set_tests_properties(prep_411_bulk_Si_gs_symmetry PROPERTIES  FIXTURES_SETUP "setup_411_bulk_Si_gs_symmetry" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/411_bulk_Si_gs_symmetry/CMakeLists.txt;1;create_test;/workspace/testsuites/411_bulk_Si_gs_symmetry/CMakeLists.txt;0;")
add_test(run_411_bulk_Si_gs_symmetry "mpiexec_411_bulk_Si_gs_symmetry.sh")
set_tests_properties(run_411_bulk_Si_gs_symmetry PROPERTIES  FIXTURES_REQUIRED "setup_411_bulk_Si_gs_symmetry" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/411_bulk_Si_gs_symmetry/CMakeLists.txt;1;create_test;/workspace/testsuites/411_bulk_Si_gs_symmetry/CMakeLists.txt;0;")
add_test(verify_411_bulk_Si_gs_symmetry "verification")
set_tests_properties(verify_411_bulk_Si_gs_symmetry PROPERTIES  FIXTURES_CLEANUP "cleanup_411_bulk_Si_gs_symmetry" FIXTURES_REQUIRED "setup_411_bulk_Si_gs_symmetry;run_411_bulk_Si_gs_symmetry" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/411_bulk_Si_gs_symmetry/CMakeLists.txt;1;create_test;/workspace/testsuites/411_bulk_Si_gs_symmetry/CMakeLists.txt;0;")
