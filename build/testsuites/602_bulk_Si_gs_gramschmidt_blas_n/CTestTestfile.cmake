# CMake generated Testfile for 
# Source directory: /workspace/testsuites/602_bulk_Si_gs_gramschmidt_blas_n
# Build directory: /workspace/build/testsuites/602_bulk_Si_gs_gramschmidt_blas_n
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_602_bulk_Si_gs_gramschmidt_blas_n "prep_602_bulk_Si_gs_gramschmidt_blas_n.sh")
set_tests_properties(prep_602_bulk_Si_gs_gramschmidt_blas_n PROPERTIES  FIXTURES_SETUP "setup_602_bulk_Si_gs_gramschmidt_blas_n" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/602_bulk_Si_gs_gramschmidt_blas_n/CMakeLists.txt;1;create_test;/workspace/testsuites/602_bulk_Si_gs_gramschmidt_blas_n/CMakeLists.txt;0;")
add_test(run_602_bulk_Si_gs_gramschmidt_blas_n "mpiexec_602_bulk_Si_gs_gramschmidt_blas_n.sh")
set_tests_properties(run_602_bulk_Si_gs_gramschmidt_blas_n PROPERTIES  FIXTURES_REQUIRED "setup_602_bulk_Si_gs_gramschmidt_blas_n" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/602_bulk_Si_gs_gramschmidt_blas_n/CMakeLists.txt;1;create_test;/workspace/testsuites/602_bulk_Si_gs_gramschmidt_blas_n/CMakeLists.txt;0;")
add_test(verify_602_bulk_Si_gs_gramschmidt_blas_n "verification")
set_tests_properties(verify_602_bulk_Si_gs_gramschmidt_blas_n PROPERTIES  FIXTURES_CLEANUP "cleanup_602_bulk_Si_gs_gramschmidt_blas_n" FIXTURES_REQUIRED "setup_602_bulk_Si_gs_gramschmidt_blas_n;run_602_bulk_Si_gs_gramschmidt_blas_n" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/602_bulk_Si_gs_gramschmidt_blas_n/CMakeLists.txt;1;create_test;/workspace/testsuites/602_bulk_Si_gs_gramschmidt_blas_n/CMakeLists.txt;0;")
