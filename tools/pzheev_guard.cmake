# User policy: PZHEEV is prohibited, including workspace-query calls.
find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(SALMON_EIGENSOLVER_GUARD "${PROJECT_SOURCE_DIR}/tools/check_forbidden_eigensolvers.py")
execute_process(
  COMMAND "${Python3_EXECUTABLE}" "${SALMON_EIGENSOLVER_GUARD}"
          --source-root "${PROJECT_SOURCE_DIR}/src"
  RESULT_VARIABLE _salmon_guard_result)
if(NOT _salmon_guard_result EQUAL 0)
  message(FATAL_ERROR "PZHEEV prohibition check failed; configure/build is blocked.")
endif()
add_custom_target(salmon_check_forbidden_eigensolvers ALL
  COMMAND "${Python3_EXECUTABLE}" "${SALMON_EIGENSOLVER_GUARD}"
          --source-root "${PROJECT_SOURCE_DIR}/src"
  VERBATIM)
function(salmon_guard_eigensolver_target target)
  add_dependencies(${target} salmon_check_forbidden_eigensolvers)
  # Only application objects, not vendor archives that export unused routines.
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND "${Python3_EXECUTABLE}" "${SALMON_EIGENSOLVER_GUARD}"
            --object-tree "${CMAKE_CURRENT_BINARY_DIR}"
            --built-executable "$<TARGET_FILE:${target}>"
    VERBATIM)
endfunction()
