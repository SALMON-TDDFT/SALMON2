# CMake generated Testfile for 
# Source directory: /workspace/testsuites/192_bulk_Si_pseudo_pspnc
# Build directory: /workspace/build/testsuites/192_bulk_Si_pseudo_pspnc
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_192_bulk_Si_pseudo_pspnc "prep_192_bulk_Si_pseudo_pspnc.sh")
set_tests_properties(prep_192_bulk_Si_pseudo_pspnc PROPERTIES  FIXTURES_SETUP "setup_192_bulk_Si_pseudo_pspnc" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/192_bulk_Si_pseudo_pspnc/CMakeLists.txt;1;create_test;/workspace/testsuites/192_bulk_Si_pseudo_pspnc/CMakeLists.txt;0;")
add_test(run_192_bulk_Si_pseudo_pspnc "mpiexec_192_bulk_Si_pseudo_pspnc.sh")
set_tests_properties(run_192_bulk_Si_pseudo_pspnc PROPERTIES  FIXTURES_REQUIRED "setup_192_bulk_Si_pseudo_pspnc" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/192_bulk_Si_pseudo_pspnc/CMakeLists.txt;1;create_test;/workspace/testsuites/192_bulk_Si_pseudo_pspnc/CMakeLists.txt;0;")
add_test(verify_192_bulk_Si_pseudo_pspnc "verification")
set_tests_properties(verify_192_bulk_Si_pseudo_pspnc PROPERTIES  FIXTURES_CLEANUP "cleanup_192_bulk_Si_pseudo_pspnc" FIXTURES_REQUIRED "setup_192_bulk_Si_pseudo_pspnc;run_192_bulk_Si_pseudo_pspnc" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/192_bulk_Si_pseudo_pspnc/CMakeLists.txt;1;create_test;/workspace/testsuites/192_bulk_Si_pseudo_pspnc/CMakeLists.txt;0;")
