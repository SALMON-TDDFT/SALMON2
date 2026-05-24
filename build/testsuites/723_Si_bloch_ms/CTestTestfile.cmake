# CMake generated Testfile for 
# Source directory: /workspace/testsuites/723_Si_bloch_ms
# Build directory: /workspace/build/testsuites/723_Si_bloch_ms
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_723_Si_bloch_ms "prep_723_Si_bloch_ms.sh")
set_tests_properties(prep_723_Si_bloch_ms PROPERTIES  FIXTURES_SETUP "setup_723_Si_bloch_ms" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/723_Si_bloch_ms/CMakeLists.txt;1;create_test;/workspace/testsuites/723_Si_bloch_ms/CMakeLists.txt;0;")
add_test(run_723_Si_bloch_ms "mpiexec_723_Si_bloch_ms.sh")
set_tests_properties(run_723_Si_bloch_ms PROPERTIES  FIXTURES_REQUIRED "setup_723_Si_bloch_ms" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/723_Si_bloch_ms/CMakeLists.txt;1;create_test;/workspace/testsuites/723_Si_bloch_ms/CMakeLists.txt;0;")
add_test(verify_723_Si_bloch_ms "verification")
set_tests_properties(verify_723_Si_bloch_ms PROPERTIES  FIXTURES_CLEANUP "cleanup_723_Si_bloch_ms" FIXTURES_REQUIRED "setup_723_Si_bloch_ms;run_723_Si_bloch_ms" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/723_Si_bloch_ms/CMakeLists.txt;1;create_test;/workspace/testsuites/723_Si_bloch_ms/CMakeLists.txt;0;")
