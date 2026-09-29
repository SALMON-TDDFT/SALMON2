set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_438_H_pbe0_mesh_rt_x PROPERTIES TIMEOUT 180 PROCESSORS 2)
set_tests_properties(run_438_H_pbe0_mesh_rt_x PROPERTIES FIXTURES_SETUP run_438_H_pbe0_mesh_rt_x)
set_tests_properties(prep_438_H_pbe0_mesh_rt_x PROPERTIES FIXTURES_REQUIRED "pbe0_dc_verified;run_424_H_pbe0_mesh_rt")
set_tests_properties(prep_438_H_pbe0_mesh_rt_x run_438_H_pbe0_mesh_rt_x verify_438_H_pbe0_mesh_rt_x
  PROPERTIES LABELS "pbe0;hybrid;block3d" ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
