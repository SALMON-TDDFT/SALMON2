add_executable(soc_wave_fft_omp_probe
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor_wave_fft.f90 probe.f90)
set(_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-457")
file(MAKE_DIRECTORY "${_modules}")
set_target_properties(soc_wave_fft_omp_probe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${_modules}")
target_include_directories(soc_wave_fft_omp_probe BEFORE PRIVATE "${_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(soc_wave_fft_omp_probe PRIVATE ${HSE_FFTW_LIBRARY})
set(SALMON_EIGENSOLVER_GUARD "${PROJECT_SOURCE_DIR}/tools/check_forbidden_eigensolvers.py")
salmon_guard_eigensolver_target(soc_wave_fft_omp_probe)
foreach(_threads 1 2 4)
 add_test(NAME 457_soc_wave_fft_omp_${_threads} COMMAND $<TARGET_FILE:soc_wave_fft_omp_probe>)
 set_tests_properties(457_soc_wave_fft_omp_${_threads} PROPERTIES TIMEOUT 30 PROCESSORS ${_threads}
  LABELS "hybrid;soc;openmp" ENVIRONMENT "OMP_NUM_THREADS=${_threads};OMP_DYNAMIC=FALSE;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(soc_wave_fft_omp_benchmark EXCLUDE_FROM_ALL
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor_wave_fft.f90 benchmark.f90)
set(_bench_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-457-bench")
file(MAKE_DIRECTORY "${_bench_modules}")
set_target_properties(soc_wave_fft_omp_benchmark PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${_bench_modules}")
target_include_directories(soc_wave_fft_omp_benchmark BEFORE PRIVATE "${_bench_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(soc_wave_fft_omp_benchmark PRIVATE ${HSE_FFTW_LIBRARY})
salmon_guard_eigensolver_target(soc_wave_fft_omp_benchmark)
add_executable(soc_wave_worker_omp_probe
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_backend.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_collinear_fft.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_exchange.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_symmetry.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor_wave_fft.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor_wave_workers.f90 worker_probe.f90)
set(_worker_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-457-worker")
file(MAKE_DIRECTORY "${_worker_modules}")
set_target_properties(soc_wave_worker_omp_probe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${_worker_modules}")
target_include_directories(soc_wave_worker_omp_probe BEFORE PRIVATE "${_worker_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(soc_wave_worker_omp_probe PRIVATE ${EXTERNAL_LIBS})
salmon_guard_eigensolver_target(soc_wave_worker_omp_probe)
foreach(_np 2 4)
 add_test(NAME 457_soc_wave_worker_omp_${_np} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${_np}
  $<TARGET_FILE:soc_wave_worker_omp_probe>)
 set_tests_properties(457_soc_wave_worker_omp_${_np} PROPERTIES TIMEOUT 60 PROCESSORS ${_np}
  LABELS "hybrid;soc;openmp;mpi" ENVIRONMENT "OMP_NUM_THREADS=4;OMP_DYNAMIC=FALSE;OPENBLAS_NUM_THREADS=1")
endforeach()
