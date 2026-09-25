create_mpi_test()
set_tests_properties(run_420_bulk_Si_hse_gs PROPERTIES FIXTURES_SETUP run_420_bulk_Si_hse_gs)
set_tests_properties(verify_420_bulk_Si_hse_gs PROPERTIES FIXTURES_SETUP hse_k4_gs_verified)
set_tests_properties(run_420_bulk_Si_hse_gs PROPERTIES TIMEOUT 1200 PROCESSORS 4 LABELS "hse;k4")

# Exercise the actual selected LAPACK, including automatically built fallbacks.
add_executable(check_hse_lapack check_lapack.f90)
target_link_libraries(check_hse_lapack PRIVATE ${EXTERNAL_LIBS})
if(EXTERNAL_PROJECT_TARGETS)
  add_dependencies(check_hse_lapack ${EXTERNAL_PROJECT_TARGETS})
endif()
add_test(NAME hse_lapack_eigenvectors COMMAND check_hse_lapack)
set_tests_properties(hse_lapack_eigenvectors PROPERTIES FIXTURES_SETUP hse_lapack_verified LABELS "hse;k4")
set_tests_properties(prep_420_bulk_Si_hse_gs PROPERTIES FIXTURES_REQUIRED hse_lapack_verified)
