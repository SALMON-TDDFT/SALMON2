set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_433_H_pbeh40_conventional_gs PROPERTIES TIMEOUT 180 PROCESSORS 2
  FIXTURES_SETUP run_433_H_pbeh40_conventional_gs)
set_tests_properties(verify_433_H_pbeh40_conventional_gs PROPERTIES FIXTURES_SETUP pbeh40_conventional_gs_verified)
set_tests_properties(prep_433_H_pbeh40_conventional_gs run_433_H_pbeh40_conventional_gs verify_433_H_pbeh40_conventional_gs
  PROPERTIES LABELS "hybrid;pbeh40;conventional;gs"
  ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
