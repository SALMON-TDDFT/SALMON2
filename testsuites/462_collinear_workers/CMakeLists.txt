add_executable(collinear_worker_probe
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_backend.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_collinear_fft.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_k_exchange.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_symmetry.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor_wave_fft.f90
  ${PROJECT_SOURCE_DIR}/src/xc/exx_spinor_wave_workers.f90 worker_probe.f90)
set(_worker_modules "${CMAKE_BINARY_DIR}/testsuite-probes/modules-462-worker")
file(MAKE_DIRECTORY "${_worker_modules}")
set_target_properties(collinear_worker_probe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/testsuite-probes"
  Fortran_MODULE_DIRECTORY "${_worker_modules}")
target_include_directories(collinear_worker_probe BEFORE PRIVATE "${_worker_modules}" "${HSE_FFTW_INCLUDE}")
target_link_libraries(collinear_worker_probe PRIVATE ${EXTERNAL_LIBS})
salmon_guard_eigensolver_target(collinear_worker_probe)
foreach(_np 2 4 8)
 add_test(NAME 462_collinear_worker_${_np} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${_np}
  $<TARGET_FILE:collinear_worker_probe>)
 set_tests_properties(462_collinear_worker_${_np} PROPERTIES TIMEOUT 60 PROCESSORS ${_np}
  LABELS "hybrid;collinear;openmp;mpi" ENVIRONMENT "OMP_NUM_THREADS=4;OMP_DYNAMIC=FALSE;OPENBLAS_NUM_THREADS=1")
endforeach()

add_test(NAME 462_collinear_native_gs_rt
 COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/verify_native.py
 --binary $<TARGET_FILE:${TARGET_NAME}>
 --output ${CMAKE_CURRENT_BINARY_DIR}/native --unique-output
 --pseudo ${PROJECT_SOURCE_DIR}/testsuites/pseudo/Si_rps.dat --threads 2)
set_tests_properties(462_collinear_native_gs_rt PROPERTIES TIMEOUT 600 PROCESSORS 16
 LABELS "hybrid;collinear;openmp;mpi;native")
