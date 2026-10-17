set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_441_H_hse_k_pair_symmetry PROPERTIES TIMEOUT 120 PROCESSORS 2)
set_tests_properties(prep_441_H_hse_k_pair_symmetry run_441_H_hse_k_pair_symmetry verify_441_H_hse_k_pair_symmetry
  PROPERTIES LABELS "hybrid;pair-symmetry"
  ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
