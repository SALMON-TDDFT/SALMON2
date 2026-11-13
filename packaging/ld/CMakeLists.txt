cmake_minimum_required(VERSION 3.20)
project(salmonpy_ld LANGUAGES Fortran)
find_package(Python REQUIRED COMPONENTS Interpreter)
get_filename_component(ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(CHECKER "${ROOT}/tools/check_forbidden_eigensolvers.py")
# Scan before compilation and on incremental builds; the source-only distribution
# may contain only the LD source, while a checkout scans all application sources.
foreach(source "${ROOT}/src" "${ROOT}/salmonpy")
  execute_process(COMMAND "${Python_EXECUTABLE}" "${CHECKER}" --source-root "${source}"
    RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "Source prohibition audit failed")
  endif()
endforeach()
add_custom_target(ld_audit ALL
  COMMAND "${Python_EXECUTABLE}" "${CHECKER}" --source-root "${ROOT}/src"
  COMMAND "${Python_EXECUTABLE}" "${CHECKER}" --source-root "${ROOT}/salmonpy")
add_library(salmonpy_ld SHARED "${ROOT}/src/maxwell/lorentz_drude.f90" "${ROOT}/salmonpy/ld_c_api.f90"
  "${ROOT}/src/ssbe/sbe_types.f90" "${ROOT}/src/ssbe/sbe_velocity_kernel.f90" "${ROOT}/salmonpy/sbe_c_api.f90")
if(APPLE)
  target_link_libraries(salmonpy_ld PRIVATE "-framework Accelerate")
else()
  set(BLA_VENDOR OpenBLAS)
  find_package(BLAS REQUIRED)
  target_link_libraries(salmonpy_ld PRIVATE BLAS::BLAS)
endif()
add_dependencies(salmonpy_ld ld_audit)
set_target_properties(salmonpy_ld PROPERTIES PREFIX "lib" WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_custom_command(TARGET salmonpy_ld POST_BUILD
  COMMAND "${Python_EXECUTABLE}" "${CHECKER}" --object-tree "${CMAKE_CURRENT_BINARY_DIR}"
  COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/hash_native.py" "$<TARGET_FILE:salmonpy_ld>")
install(TARGETS salmonpy_ld LIBRARY DESTINATION salmonpy/_native RUNTIME DESTINATION salmonpy/_native)
install(FILES "$<TARGET_FILE:salmonpy_ld>.sha256" DESTINATION salmonpy/_native)
install(FILES "${ROOT}/salmonpy/__init__.py" "${ROOT}/salmonpy/ld.py"
  "${ROOT}/salmonpy/fields.py" "${ROOT}/salmonpy/geometry.py"
  "${ROOT}/salmonpy/analysis.py" "${ROOT}/salmonpy/fdtd.py" "${ROOT}/salmonpy/crystal.py" "${ROOT}/salmonpy/sbe.py"
  "${ROOT}/salmonpy/WORKFLOWS.md" "${ROOT}/salmonpy/PACKAGING.md" "${ROOT}/salmonpy/README.md" "${ROOT}/salmonpy/AI_USAGE.md" "${ROOT}/salmonpy/INPUT_CONTRACT.md"
  DESTINATION salmonpy)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/licenses/" DESTINATION salmonpy/licenses)
