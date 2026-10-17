set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_440_H_pbeh40_symmetry_gs PROPERTIES TIMEOUT 120 PROCESSORS 1)
set_tests_properties(prep_440_H_pbeh40_symmetry_gs run_440_H_pbeh40_symmetry_gs verify_440_H_pbeh40_symmetry_gs
  PROPERTIES LABELS "hybrid;symmetry")
