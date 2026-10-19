add_test(NAME 446_checkpoint_elapsed COMMAND python3
  ${CMAKE_CURRENT_SOURCE_DIR}/verification
  ${CMAKE_Fortran_COMPILER} ${PROJECT_SOURCE_DIR}/src/rt/main_tddft.f90)
set_tests_properties(446_checkpoint_elapsed PROPERTIES LABELS "hybrid;checkpoint" TIMEOUT 30)
