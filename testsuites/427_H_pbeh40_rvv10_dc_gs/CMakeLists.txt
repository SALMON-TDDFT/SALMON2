set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_427_H_pbeh40_rvv10_dc_gs PROPERTIES TIMEOUT 180 PROCESSORS 2
  FIXTURES_SETUP run_427_H_pbeh40_rvv10_dc_gs)
set_tests_properties(verify_427_H_pbeh40_rvv10_dc_gs PROPERTIES FIXTURES_SETUP pbeh40_rvv10_dc_verified)
set_tests_properties(prep_427_H_pbeh40_rvv10_dc_gs run_427_H_pbeh40_rvv10_dc_gs verify_427_H_pbeh40_rvv10_dc_gs
  PROPERTIES LABELS "hybrid;pbeh40_rvv10;gs"
  ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
