set(TEST_COMMAND "${CMAKE_COMMAND} -E env FC=${CMAKE_Fortran_COMPILER} ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_periodic_grid.py")
create_test()
set_tests_properties(run_452_symmetry_periodic_grid PROPERTIES FIXTURES_SETUP run_452_symmetry_periodic_grid TIMEOUT 90 PROCESSORS 2)
set_tests_properties(prep_452_symmetry_periodic_grid run_452_symmetry_periodic_grid verify_452_symmetry_periodic_grid PROPERTIES LABELS "symmetry;periodic;unit")
