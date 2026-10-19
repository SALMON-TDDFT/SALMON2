# Uppercase is normalized; enabling spinor tiles without SOC must fail early.
set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_444_spinor_pair_input PROPERTIES
  FIXTURES_SETUP run_444_spinor_pair_input TIMEOUT 60 PROCESSORS 2)
set_tests_properties(prep_444_spinor_pair_input run_444_spinor_pair_input verify_444_spinor_pair_input
  PROPERTIES LABELS "hybrid;spinor-pair" ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
