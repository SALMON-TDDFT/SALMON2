set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_424_H_pbe0_mesh_rt PROPERTIES TIMEOUT 180 PROCESSORS 2)
set_tests_properties(run_424_H_pbe0_mesh_rt PROPERTIES FIXTURES_SETUP run_424_H_pbe0_mesh_rt)
set_tests_properties(prep_424_H_pbe0_mesh_rt PROPERTIES FIXTURES_REQUIRED pbe0_dc_verified)
set_tests_properties(prep_424_H_pbe0_mesh_rt run_424_H_pbe0_mesh_rt verify_424_H_pbe0_mesh_rt
  PROPERTIES LABELS "pbe0;hybrid")
