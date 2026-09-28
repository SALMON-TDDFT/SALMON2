create_mpi_test()
if(USE_MPI)
  set_tests_properties(run_422_H_dcdft_hse PROPERTIES FIXTURES_SETUP run_422_H_dcdft_hse)
endif()
