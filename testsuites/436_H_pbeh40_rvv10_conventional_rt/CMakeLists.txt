set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_436_H_pbeh40_rvv10_conventional_rt PROPERTIES TIMEOUT 180 PROCESSORS 2
  FIXTURES_SETUP run_436_H_pbeh40_rvv10_conventional_rt)
set_tests_properties(prep_436_H_pbeh40_rvv10_conventional_rt PROPERTIES FIXTURES_REQUIRED pbeh40_rvv10_conventional_gs_verified)
set_tests_properties(prep_436_H_pbeh40_rvv10_conventional_rt run_436_H_pbeh40_rvv10_conventional_rt verify_436_H_pbeh40_rvv10_conventional_rt
  PROPERTIES LABELS "hybrid;pbeh40_rvv10;conventional;rt"
  ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
