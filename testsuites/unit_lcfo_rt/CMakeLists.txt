# Separate module directory prevents a parallel build from racing the main target.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
add_executable(check_lcfo_rt_core
  "${CMAKE_SOURCE_DIR}/src/rt/lcfo_rt_core.f90" probe.f90)
set_target_properties(check_lcfo_rt_core PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(check_lcfo_rt_core PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_link_libraries(check_lcfo_rt_core PRIVATE ${EXTERNAL_LIBS})
if(EXTERNAL_PROJECT_TARGETS)
  add_dependencies(check_lcfo_rt_core ${EXTERNAL_PROJECT_TARGETS})
endif()
add_test(NAME lcfo_rt_core COMMAND check_lcfo_rt_core)
set_tests_properties(lcfo_rt_core PROPERTIES LABELS "lcfo;rt;unit"
  ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1" TIMEOUT 30)
