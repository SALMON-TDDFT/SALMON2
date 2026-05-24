# CMake generated Testfile for 
# Source directory: /workspace/testsuites/413_bulk_Si_rt_response_symmetry_yyn
# Build directory: /workspace/build/testsuites/413_bulk_Si_rt_response_symmetry_yyn
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_413_bulk_Si_rt_response_symmetry_yyn "prep_413_bulk_Si_rt_response_symmetry_yyn.sh")
set_tests_properties(prep_413_bulk_Si_rt_response_symmetry_yyn PROPERTIES  FIXTURES_SETUP "setup_413_bulk_Si_rt_response_symmetry_yyn" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/413_bulk_Si_rt_response_symmetry_yyn/CMakeLists.txt;1;create_test;/workspace/testsuites/413_bulk_Si_rt_response_symmetry_yyn/CMakeLists.txt;0;")
add_test(run_413_bulk_Si_rt_response_symmetry_yyn "mpiexec_413_bulk_Si_rt_response_symmetry_yyn.sh")
set_tests_properties(run_413_bulk_Si_rt_response_symmetry_yyn PROPERTIES  FIXTURES_REQUIRED "setup_413_bulk_Si_rt_response_symmetry_yyn" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/413_bulk_Si_rt_response_symmetry_yyn/CMakeLists.txt;1;create_test;/workspace/testsuites/413_bulk_Si_rt_response_symmetry_yyn/CMakeLists.txt;0;")
add_test(verify_413_bulk_Si_rt_response_symmetry_yyn "verification")
set_tests_properties(verify_413_bulk_Si_rt_response_symmetry_yyn PROPERTIES  FIXTURES_CLEANUP "cleanup_413_bulk_Si_rt_response_symmetry_yyn" FIXTURES_REQUIRED "setup_413_bulk_Si_rt_response_symmetry_yyn;run_413_bulk_Si_rt_response_symmetry_yyn" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/413_bulk_Si_rt_response_symmetry_yyn/CMakeLists.txt;1;create_test;/workspace/testsuites/413_bulk_Si_rt_response_symmetry_yyn/CMakeLists.txt;0;")
