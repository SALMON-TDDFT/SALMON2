add_executable(fft_consensus_probe
  ${PROJECT_SOURCE_DIR}/src/misc/nvtx_wrapper.f90
  ${PROJECT_SOURCE_DIR}/src/parallel/communication.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_backend.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_collinear_fft.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_exchange.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_fft_consensus.f90
  probe.f90)
set(_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-fft-consensus")
file(MAKE_DIRECTORY "${_modules}")
set_target_properties(fft_consensus_probe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${_modules}")
target_include_directories(fft_consensus_probe BEFORE PRIVATE "${_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(fft_consensus_probe PRIVATE ${EXTERNAL_LIBS})
set(SALMON_EIGENSOLVER_GUARD "${PROJECT_SOURCE_DIR}/tools/check_forbidden_eigensolvers.py")
salmon_guard_eigensolver_target(fft_consensus_probe)
foreach(_np 2 4 8)
  add_test(NAME fft_consensus_${_np} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${_np} $<TARGET_FILE:fft_consensus_probe>)
  set_tests_properties(fft_consensus_${_np} PROPERTIES TIMEOUT 90 PROCESSORS ${_np} LABELS "hybrid;mpi;fft;unit"
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
endforeach()
