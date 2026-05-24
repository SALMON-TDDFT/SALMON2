# CMake generated Testfile for 
# Source directory: /workspace/testsuites/712_bulk_Si_nonorthogonal_lattice_pre
# Build directory: /workspace/build/testsuites/712_bulk_Si_nonorthogonal_lattice_pre
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_712_bulk_Si_nonorthogonal_lattice_pre "prep_712_bulk_Si_nonorthogonal_lattice_pre.sh")
set_tests_properties(prep_712_bulk_Si_nonorthogonal_lattice_pre PROPERTIES  FIXTURES_SETUP "setup_712_bulk_Si_nonorthogonal_lattice_pre" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/712_bulk_Si_nonorthogonal_lattice_pre/CMakeLists.txt;1;create_test;/workspace/testsuites/712_bulk_Si_nonorthogonal_lattice_pre/CMakeLists.txt;0;")
add_test(run_712_bulk_Si_nonorthogonal_lattice_pre "mpiexec_712_bulk_Si_nonorthogonal_lattice_pre.sh")
set_tests_properties(run_712_bulk_Si_nonorthogonal_lattice_pre PROPERTIES  FIXTURES_REQUIRED "setup_712_bulk_Si_nonorthogonal_lattice_pre" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/712_bulk_Si_nonorthogonal_lattice_pre/CMakeLists.txt;1;create_test;/workspace/testsuites/712_bulk_Si_nonorthogonal_lattice_pre/CMakeLists.txt;0;")
add_test(verify_712_bulk_Si_nonorthogonal_lattice_pre "verification")
set_tests_properties(verify_712_bulk_Si_nonorthogonal_lattice_pre PROPERTIES  FIXTURES_CLEANUP "cleanup_712_bulk_Si_nonorthogonal_lattice_pre" FIXTURES_REQUIRED "setup_712_bulk_Si_nonorthogonal_lattice_pre;run_712_bulk_Si_nonorthogonal_lattice_pre" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/712_bulk_Si_nonorthogonal_lattice_pre/CMakeLists.txt;1;create_test;/workspace/testsuites/712_bulk_Si_nonorthogonal_lattice_pre/CMakeLists.txt;0;")
