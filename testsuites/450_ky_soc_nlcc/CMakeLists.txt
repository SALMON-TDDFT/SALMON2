set(TEST_COMMAND "${CMAKE_COMMAND} -E env FC=${CMAKE_Fortran_COMPILER} ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.py")
create_test()
set_tests_properties(run_450_ky_soc_nlcc PROPERTIES FIXTURES_SETUP run_450_ky_soc_nlcc TIMEOUT 180 PROCESSORS 1)
set_tests_properties(prep_450_ky_soc_nlcc run_450_ky_soc_nlcc verify_450_ky_soc_nlcc PROPERTIES LABELS "pseudo;soc;nlcc;generator;unit")
