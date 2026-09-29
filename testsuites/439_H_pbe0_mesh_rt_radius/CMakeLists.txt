set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
create_mpi_test()
set_tests_properties(run_439_H_pbe0_mesh_rt_radius PROPERTIES TIMEOUT 180 PROCESSORS 2)
set_tests_properties(run_439_H_pbe0_mesh_rt_radius PROPERTIES FIXTURES_SETUP run_439_H_pbe0_mesh_rt_radius)
set_tests_properties(prep_439_H_pbe0_mesh_rt_radius PROPERTIES FIXTURES_REQUIRED pbe0_dc_verified)
set_tests_properties(prep_439_H_pbe0_mesh_rt_radius run_439_H_pbe0_mesh_rt_radius verify_439_H_pbe0_mesh_rt_radius
  PROPERTIES LABELS "pbe0;hybrid;block3d" ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
