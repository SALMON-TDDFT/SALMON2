# Exercise the actual selected LAPACK, including automatically built fallbacks.
add_executable(check_hse_lapack check_lapack.f90)
target_link_libraries(check_hse_lapack PRIVATE ${EXTERNAL_LIBS})
if(EXTERNAL_PROJECT_TARGETS)
  add_dependencies(check_hse_lapack ${EXTERNAL_PROJECT_TARGETS})
endif()
add_test(NAME hse_lapack_eigenvectors COMMAND check_hse_lapack)
set_tests_properties(hse_lapack_eigenvectors PROPERTIES FIXTURES_SETUP hse_lapack_verified LABELS "hse;unit")
