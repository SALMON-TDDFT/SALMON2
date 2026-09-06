# The policy has no MPI dependency; exercise large rank counts without launching them.
add_executable(test_automatic_distribution
  ${CMAKE_SOURCE_DIR}/src/parallel/automatic_distribution.f90
  test_automatic_distribution.f90)
set_target_properties(test_automatic_distribution PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(test_automatic_distribution PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME automatic_distribution COMMAND test_automatic_distribution)
