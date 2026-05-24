# CMake generated Testfile for 
# Source directory: /workspace/testsuites/195_bulk_Si_pseudo_upf
# Build directory: /workspace/build/testsuites/195_bulk_Si_pseudo_upf
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_195_bulk_Si_pseudo_upf "prep_195_bulk_Si_pseudo_upf.sh")
set_tests_properties(prep_195_bulk_Si_pseudo_upf PROPERTIES  FIXTURES_SETUP "setup_195_bulk_Si_pseudo_upf" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/195_bulk_Si_pseudo_upf/CMakeLists.txt;1;create_test;/workspace/testsuites/195_bulk_Si_pseudo_upf/CMakeLists.txt;0;")
add_test(run_195_bulk_Si_pseudo_upf "mpiexec_195_bulk_Si_pseudo_upf.sh")
set_tests_properties(run_195_bulk_Si_pseudo_upf PROPERTIES  FIXTURES_REQUIRED "setup_195_bulk_Si_pseudo_upf" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/195_bulk_Si_pseudo_upf/CMakeLists.txt;1;create_test;/workspace/testsuites/195_bulk_Si_pseudo_upf/CMakeLists.txt;0;")
add_test(verify_195_bulk_Si_pseudo_upf "verification")
set_tests_properties(verify_195_bulk_Si_pseudo_upf PROPERTIES  FIXTURES_CLEANUP "cleanup_195_bulk_Si_pseudo_upf" FIXTURES_REQUIRED "setup_195_bulk_Si_pseudo_upf;run_195_bulk_Si_pseudo_upf" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/195_bulk_Si_pseudo_upf/CMakeLists.txt;1;create_test;/workspace/testsuites/195_bulk_Si_pseudo_upf/CMakeLists.txt;0;")
