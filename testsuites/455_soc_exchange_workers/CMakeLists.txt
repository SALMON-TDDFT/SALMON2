set(TEST_COMMAND "python3 ${CMAKE_CURRENT_SOURCE_DIR}/verify_native.py --binary ${CMAKE_BINARY_DIR}/${TARGET_NAME} --output ${CMAKE_CURRENT_BINARY_DIR}/comparison --mpiexec ${MPIEXEC}")
create_mpi_test()
set_tests_properties(run_455_soc_exchange_workers PROPERTIES FIXTURES_SETUP run_455_soc_exchange_workers TIMEOUT 600)
set_tests_properties(prep_455_soc_exchange_workers run_455_soc_exchange_workers verify_455_soc_exchange_workers
  PROPERTIES LABELS "hybrid;soc;exchange-workers" PROCESSORS 12)
