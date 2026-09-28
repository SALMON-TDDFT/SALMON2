set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_437_H_pbe0_conventional_pulse PROPERTIES TIMEOUT 180 PROCESSORS 2
  FIXTURES_SETUP run_437_H_pbe0_conventional_pulse)
set_tests_properties(prep_437_H_pbe0_conventional_pulse PROPERTIES FIXTURES_REQUIRED pbe0_conventional_gs_verified)
set_tests_properties(prep_437_H_pbe0_conventional_pulse run_437_H_pbe0_conventional_pulse verify_437_H_pbe0_conventional_pulse
  PROPERTIES LABELS "hybrid;pbe0;conventional;rt"
  ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
