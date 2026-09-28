set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_429_H_pbeh40_mesh_ehrenfest PROPERTIES TIMEOUT 180 PROCESSORS 2
  FIXTURES_SETUP run_429_H_pbeh40_mesh_ehrenfest)
set_tests_properties(prep_429_H_pbeh40_mesh_ehrenfest PROPERTIES FIXTURES_REQUIRED pbeh40_dc_verified)
set_tests_properties(prep_429_H_pbeh40_mesh_ehrenfest run_429_H_pbeh40_mesh_ehrenfest verify_429_H_pbeh40_mesh_ehrenfest
  PROPERTIES LABELS "hybrid;pbeh40;rt;md"
  ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
