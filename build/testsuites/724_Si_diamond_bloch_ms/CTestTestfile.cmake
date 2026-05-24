# CMake generated Testfile for 
# Source directory: /workspace/testsuites/724_Si_diamond_bloch_ms
# Build directory: /workspace/build/testsuites/724_Si_diamond_bloch_ms
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_724_Si_diamond_bloch_ms "prep_724_Si_diamond_bloch_ms.sh")
set_tests_properties(prep_724_Si_diamond_bloch_ms PROPERTIES  FIXTURES_SETUP "setup_724_Si_diamond_bloch_ms" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/724_Si_diamond_bloch_ms/CMakeLists.txt;1;create_test;/workspace/testsuites/724_Si_diamond_bloch_ms/CMakeLists.txt;0;")
add_test(run_724_Si_diamond_bloch_ms "mpiexec_724_Si_diamond_bloch_ms.sh")
set_tests_properties(run_724_Si_diamond_bloch_ms PROPERTIES  FIXTURES_REQUIRED "setup_724_Si_diamond_bloch_ms" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/724_Si_diamond_bloch_ms/CMakeLists.txt;1;create_test;/workspace/testsuites/724_Si_diamond_bloch_ms/CMakeLists.txt;0;")
add_test(verify_724_Si_diamond_bloch_ms "verification")
set_tests_properties(verify_724_Si_diamond_bloch_ms PROPERTIES  FIXTURES_CLEANUP "cleanup_724_Si_diamond_bloch_ms" FIXTURES_REQUIRED "setup_724_Si_diamond_bloch_ms;run_724_Si_diamond_bloch_ms" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/724_Si_diamond_bloch_ms/CMakeLists.txt;1;create_test;/workspace/testsuites/724_Si_diamond_bloch_ms/CMakeLists.txt;0;")
