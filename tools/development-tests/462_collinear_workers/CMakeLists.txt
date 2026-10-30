add_test(NAME 462_worker_layout_oracle COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/test_layout.py
  --binary $<TARGET_FILE:worker_layout_probe>)
set_tests_properties(462_worker_layout_oracle PROPERTIES TIMEOUT 60 LABELS "hybrid;layout")

add_test(NAME 462_collinear_native_gs_rt
 COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/verify_native.py
 --binary $<TARGET_FILE:${TARGET_NAME}>
 --output ${CMAKE_CURRENT_BINARY_DIR}/native --unique-output
 --pseudo ${PROJECT_SOURCE_DIR}/testsuites/pseudo/Si_rps.dat --threads 2)
set_tests_properties(462_collinear_native_gs_rt PROPERTIES TIMEOUT 600 PROCESSORS 16
 LABELS "hybrid;collinear;openmp;mpi;native")

add_test(NAME 462_balanced_strip_oracle COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/test_balanced.py
 --binary $<TARGET_FILE:balanced_strip_probe>)
set_tests_properties(462_balanced_strip_oracle PROPERTIES TIMEOUT 60 LABELS "hybrid;layout;balanced")

# Reuse identical GS/RT numerical and all-rank checks for both no-PC modes.
foreach(_factor off learned)
 add_test(NAME 462_collinear_native_no_pc_${_factor}
  COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/verify_native.py
  --binary $<TARGET_FILE:${TARGET_NAME}>
  --output ${CMAKE_CURRENT_BINARY_DIR}/native-no-pc-${_factor} --unique-output
  --pseudo ${PROJECT_SOURCE_DIR}/testsuites/pseudo/Si_rps.dat --threads 2
  --balanced --no-pc --factor-mode ${_factor})
 set_tests_properties(462_collinear_native_no_pc_${_factor} PROPERTIES TIMEOUT 600 PROCESSORS 82
  LABELS "hybrid;collinear;openmp;mpi;native;balanced;no-pc")
endforeach()
