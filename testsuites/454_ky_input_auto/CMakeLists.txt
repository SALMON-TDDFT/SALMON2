set(TEST_COMMAND "${CMAKE_COMMAND} -E env FC=${CMAKE_Fortran_COMPILER} ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_auto.py")
create_test()
set_tests_properties(run_454_ky_input_auto PROPERTIES FIXTURES_SETUP run_454_ky_input_auto TIMEOUT 90 PROCESSORS 1)
set_tests_properties(prep_454_ky_input_auto run_454_ky_input_auto verify_454_ky_input_auto PROPERTIES LABELS "pseudo;soc;input;unit")
if(USE_MPI)
  add_test(NAME ky_auto_input_smoke COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../451_ky_reader_smoke/run_smoke.py
    --binary ${CMAKE_BINARY_DIR}/${TARGET_NAME} --auto-input --mpi ${MPIEXEC_EXECUTABLE} --np-flag=${MPIEXEC_NUMPROC_FLAG})
else()
  add_test(NAME ky_auto_input_smoke COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../451_ky_reader_smoke/run_smoke.py
    --binary ${CMAKE_BINARY_DIR}/${TARGET_NAME} --auto-input)
endif()
set_tests_properties(ky_auto_input_smoke PROPERTIES LABELS "pseudo;soc;input;integration" TIMEOUT 480 PROCESSORS 2)
