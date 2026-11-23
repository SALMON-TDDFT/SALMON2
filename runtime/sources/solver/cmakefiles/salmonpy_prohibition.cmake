find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(SALMONPY_CHECKER "${CMAKE_SOURCE_DIR}/tools/check_forbidden_eigensolvers.py")
execute_process(COMMAND "${Python3_EXECUTABLE}" "${SALMONPY_CHECKER}" --source-root "${CMAKE_SOURCE_DIR}/src"
  RESULT_VARIABLE SALMONPY_AUDIT_STATUS)
if(NOT SALMONPY_AUDIT_STATUS EQUAL 0)
  message(FATAL_ERROR "Required eigensolver source audit failed")
endif()
add_custom_target(salmonpy_prohibition_sources ALL
  COMMAND "${Python3_EXECUTABLE}" "${SALMONPY_CHECKER}" --source-root "${CMAKE_SOURCE_DIR}/src" VERBATIM)
