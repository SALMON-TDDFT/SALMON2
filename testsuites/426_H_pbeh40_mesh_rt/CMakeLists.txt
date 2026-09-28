set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_426_H_pbeh40_mesh_rt PROPERTIES TIMEOUT 180 PROCESSORS 2
  FIXTURES_SETUP run_426_H_pbeh40_mesh_rt)
set_tests_properties(prep_426_H_pbeh40_mesh_rt PROPERTIES FIXTURES_REQUIRED pbeh40_dc_verified)
set_tests_properties(prep_426_H_pbeh40_mesh_rt run_426_H_pbeh40_mesh_rt verify_426_H_pbeh40_mesh_rt
  PROPERTIES LABELS "hybrid;pbeh40;rt"
  ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
