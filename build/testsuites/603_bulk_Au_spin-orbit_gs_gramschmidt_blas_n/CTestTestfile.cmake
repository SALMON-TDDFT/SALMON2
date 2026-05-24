# CMake generated Testfile for 
# Source directory: /workspace/testsuites/603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n
# Build directory: /workspace/build/testsuites/603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(prep_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n "prep_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n.sh")
set_tests_properties(prep_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n PROPERTIES  FIXTURES_SETUP "setup_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;32;add_test;/workspace/testsuites/603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n/CMakeLists.txt;1;create_test;/workspace/testsuites/603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n/CMakeLists.txt;0;")
add_test(run_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n "mpiexec_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n.sh")
set_tests_properties(run_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n PROPERTIES  FIXTURES_REQUIRED "setup_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n" RUN_SERIAL "ENABLE" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;33;add_test;/workspace/testsuites/603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n/CMakeLists.txt;1;create_test;/workspace/testsuites/603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n/CMakeLists.txt;0;")
add_test(verify_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n "verification")
set_tests_properties(verify_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n PROPERTIES  FIXTURES_CLEANUP "cleanup_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n" FIXTURES_REQUIRED "setup_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n;run_603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n" _BACKTRACE_TRIPLES "/workspace/cmakefiles/create_test.cmake;34;add_test;/workspace/testsuites/603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n/CMakeLists.txt;1;create_test;/workspace/testsuites/603_bulk_Au_spin-orbit_gs_gramschmidt_blas_n/CMakeLists.txt;0;")
