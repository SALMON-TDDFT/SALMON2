# CMake generated Testfile for 
# Source directory: /workspace/testsuites/725_Si_bloch_3rd_order
# Build directory: /workspace/build/testsuites/725_Si_bloch_3rd_order
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_725_Si_bloch_3rd_order "prep_725_Si_bloch_3rd_order.sh")
set_tests_properties(prep_725_Si_bloch_3rd_order PROPERTIES  FIXTURES_SETUP "setup_725_Si_bloch_3rd_order" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/725_Si_bloch_3rd_order/CMakeLists.txt;1;create_test;/workspace/testsuites/725_Si_bloch_3rd_order/CMakeLists.txt;0;")
add_test(run_725_Si_bloch_3rd_order "mpiexec_725_Si_bloch_3rd_order.sh")
set_tests_properties(run_725_Si_bloch_3rd_order PROPERTIES  FIXTURES_REQUIRED "setup_725_Si_bloch_3rd_order" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/725_Si_bloch_3rd_order/CMakeLists.txt;1;create_test;/workspace/testsuites/725_Si_bloch_3rd_order/CMakeLists.txt;0;")
add_test(verify_725_Si_bloch_3rd_order "verification")
set_tests_properties(verify_725_Si_bloch_3rd_order PROPERTIES  FIXTURES_CLEANUP "cleanup_725_Si_bloch_3rd_order" FIXTURES_REQUIRED "setup_725_Si_bloch_3rd_order;run_725_Si_bloch_3rd_order" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/725_Si_bloch_3rd_order/CMakeLists.txt;1;create_test;/workspace/testsuites/725_Si_bloch_3rd_order/CMakeLists.txt;0;")
