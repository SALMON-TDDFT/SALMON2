set(SALMON_EIGENSOLVER_GUARD "${PROJECT_SOURCE_DIR}/tools/check_forbidden_eigensolvers.py")
add_executable(collinear_fft_omp_probe ${PROJECT_SOURCE_DIR}/src/xc/exx_collinear_fft.f90 probe.f90)
set(_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-458")
file(MAKE_DIRECTORY "${_modules}")
set_target_properties(collinear_fft_omp_probe PROPERTIES
 RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes" Fortran_MODULE_DIRECTORY "${_modules}")
target_include_directories(collinear_fft_omp_probe BEFORE PRIVATE "${_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(collinear_fft_omp_probe PRIVATE ${HSE_FFTW_LIBRARY})
salmon_guard_eigensolver_target(collinear_fft_omp_probe)
foreach(_threads 1 2 4)
 add_test(NAME 458_collinear_fft_omp_${_threads} COMMAND $<TARGET_FILE:collinear_fft_omp_probe>)
 set_tests_properties(458_collinear_fft_omp_${_threads} PROPERTIES TIMEOUT 30 PROCESSORS ${_threads}
 LABELS "hybrid;collinear;openmp" ENVIRONMENT "OMP_NUM_THREADS=${_threads};OMP_DYNAMIC=FALSE;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(collinear_action_omp_probe
 ${PROJECT_SOURCE_DIR}/src/misc/nvtx_wrapper.f90
 ${PROJECT_SOURCE_DIR}/src/parallel/communication.f90
 ${PROJECT_SOURCE_DIR}/src/xc/exx_k_backend.f90
 ${PROJECT_SOURCE_DIR}/src/xc/exx_collinear_fft.f90
 ${PROJECT_SOURCE_DIR}/src/xc/exx_k_exchange.f90
 ${PROJECT_SOURCE_DIR}/src/xc/exx_symmetry.f90
 ${PROJECT_SOURCE_DIR}/src/xc/exx_k_pairs.f90
 ${PROJECT_SOURCE_DIR}/src/xc/exx_k_orbits.f90 action_probe.f90)
set(_action_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-458-action")
file(MAKE_DIRECTORY "${_action_modules}")
set_target_properties(collinear_action_omp_probe PROPERTIES
 RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes" Fortran_MODULE_DIRECTORY "${_action_modules}")
target_include_directories(collinear_action_omp_probe BEFORE PRIVATE "${_action_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(collinear_action_omp_probe PRIVATE ${EXTERNAL_LIBS})
salmon_guard_eigensolver_target(collinear_action_omp_probe)
foreach(_np 2 4 12)
 add_test(NAME 458_collinear_action_omp_${_np} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${_np} $<TARGET_FILE:collinear_action_omp_probe>)
 set_tests_properties(458_collinear_action_omp_${_np} PROPERTIES TIMEOUT 90 PROCESSORS ${_np}
 LABELS "hybrid;collinear;openmp;mpi" ENVIRONMENT "OMP_NUM_THREADS=4;OMP_DYNAMIC=FALSE;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(collinear_fft_omp_benchmark EXCLUDE_FROM_ALL
 ${PROJECT_SOURCE_DIR}/src/xc/exx_collinear_fft.f90 benchmark.f90)
set(_bench_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-458-bench")
file(MAKE_DIRECTORY "${_bench_modules}")
set_target_properties(collinear_fft_omp_benchmark PROPERTIES
 RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes" Fortran_MODULE_DIRECTORY "${_bench_modules}")
target_include_directories(collinear_fft_omp_benchmark BEFORE PRIVATE "${_bench_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(collinear_fft_omp_benchmark PRIVATE ${HSE_FFTW_LIBRARY})
salmon_guard_eigensolver_target(collinear_fft_omp_benchmark)
