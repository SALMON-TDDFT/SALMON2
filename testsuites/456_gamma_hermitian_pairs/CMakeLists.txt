add_executable(gamma_hermitian_pairs_probe
  ${PROJECT_SOURCE_DIR}/src/misc/nvtx_wrapper.f90
  ${PROJECT_SOURCE_DIR}/src/parallel/communication.f90
  ${PROJECT_SOURCE_DIR}/src/xc/fftw_blocks.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_gamma_pairs.f90
  probe.f90)
set(_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-456")
file(MAKE_DIRECTORY "${_modules}")
set_target_properties(gamma_hermitian_pairs_probe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${_modules}")
target_include_directories(gamma_hermitian_pairs_probe BEFORE PRIVATE "${_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(gamma_hermitian_pairs_probe PRIVATE ${EXTERNAL_LIBS})
set(SALMON_EIGENSOLVER_GUARD "${PROJECT_SOURCE_DIR}/tools/check_forbidden_eigensolvers.py")
salmon_guard_eigensolver_target(gamma_hermitian_pairs_probe)
foreach(_np 1 2 4 8)
  add_test(NAME 456_gamma_pairs_${_np} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${_np}
    $<TARGET_FILE:gamma_hermitian_pairs_probe>)
  set_tests_properties(456_gamma_pairs_${_np} PROPERTIES TIMEOUT 90 PROCESSORS ${_np}
    LABELS "hybrid;dc;hermitian" ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
foreach(_mode spatial nan negative inf weight_mismatch source_nan source_inf)
  add_test(NAME 456_gamma_pairs_${_mode} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
    $<TARGET_FILE:gamma_hermitian_pairs_probe> ${_mode})
  set_tests_properties(456_gamma_pairs_${_mode} PROPERTIES TIMEOUT 90 PROCESSORS 4
    LABELS "hybrid;dc;hermitian" ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
