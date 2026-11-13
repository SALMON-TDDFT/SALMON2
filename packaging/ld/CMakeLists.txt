cmake_minimum_required(VERSION 3.20)
project(salmonpy_ld LANGUAGES Fortran)
find_package(Python REQUIRED COMPONENTS Interpreter)
get_filename_component(ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(CHECKER "${ROOT}/tools/check_forbidden_eigensolvers.py")
# Audit the bundled kernels; the optional runtime is audited after preparation.
foreach(source "${ROOT}/native" "${ROOT}/salmonpy")
  execute_process(COMMAND "${Python_EXECUTABLE}" "${CHECKER}" --source-root "${source}"
    RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "Source prohibition audit failed")
  endif()
endforeach()
add_custom_target(ld_audit ALL
  COMMAND "${Python_EXECUTABLE}" "${CHECKER}" --source-root "${ROOT}/native"
  COMMAND "${Python_EXECUTABLE}" "${CHECKER}" --source-root "${ROOT}/salmonpy")
add_library(salmonpy_ld SHARED "${ROOT}/native/lorentz_drude.f90" "${ROOT}/salmonpy/ld_c_api.f90"
  "${ROOT}/native/sbe_types.f90" "${ROOT}/native/sbe_velocity_kernel.f90" "${ROOT}/salmonpy/sbe_c_api.f90")
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
  "${ROOT}/salmonpy/salmon_vr.py" "${ROOT}/salmonpy/SALMON_VR.md" "${ROOT}/salmonpy/SALMON_VR.en.md" "${ROOT}/salmonpy/doctor.py" "${ROOT}/salmonpy/AGENTS.md" "${ROOT}/salmonpy/manual.py" "${ROOT}/salmonpy/manual.json" "${ROOT}/salmonpy/ase_calculator.py" "${ROOT}/salmonpy/dft.py" "${ROOT}/salmonpy/ase.py" "${ROOT}/salmonpy/bands.py" "${ROOT}/salmonpy/analysis.py" "${ROOT}/salmonpy/outputs.py" "${ROOT}/salmonpy/fdtd.py" "${ROOT}/salmonpy/crystal.py" "${ROOT}/salmonpy/sbe.py" "${ROOT}/salmonpy/tddft.py"
  "${ROOT}/salmonpy/ASE_CALCULATOR.md" "${ROOT}/salmonpy/DFT.md" "${ROOT}/salmonpy/ASE.md" "${ROOT}/salmonpy/BANDS.md" "${ROOT}/salmonpy/WORKFLOWS.md" "${ROOT}/salmonpy/OUTPUTS.md" "${ROOT}/salmonpy/PACKAGING.md" "${ROOT}/salmonpy/README.md" "${ROOT}/salmonpy/README.en.md" "${ROOT}/salmonpy/AI_USAGE.md" "${ROOT}/salmonpy/INPUT_CONTRACT.md"
  DESTINATION salmonpy)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/licenses/" DESTINATION salmonpy/licenses)

install(DIRECTORY "${ROOT}/salmonpy/examples/" DESTINATION salmonpy/examples
  FILES_MATCHING PATTERN "*.py" PATTERN "*.md" PATTERN "__pycache__" EXCLUDE)

install(FILES "${ROOT}/NOTICE" "${ROOT}/LICENSE" "${ROOT}/LICENSE.THIRD-PARTY" DESTINATION salmonpy/licenses)

# Explicit runtime setup must also work from an installed wheel, without Git history.
install(FILES "${ROOT}/salmonpy/runtime.py" "${ROOT}/salmonpy/runtime_catalog.json"
  "${ROOT}/salmonpy/FDTD_MATERIALS.md" "${ROOT}/salmonpy/FDTD_MATERIALS.en.md"
  "${ROOT}/salmonpy/RUNTIME.md" "${ROOT}/salmonpy/RUNTIME.en.md" DESTINATION salmonpy)
install(FILES "${ROOT}/tools/prepare_salmon_runtime.py" "${ROOT}/tools/check_forbidden_eigensolvers.py"
  DESTINATION salmonpy/_runtime_source/tools)
install(DIRECTORY "${ROOT}/runtime/" DESTINATION salmonpy/_runtime_source/runtime
  FILES_MATCHING PATTERN "*.json" PATTERN "*.patch")
install(DIRECTORY "${ROOT}/native/" DESTINATION salmonpy/_runtime_source/native
  FILES_MATCHING PATTERN "*.f90")

install(FILES "${ROOT}/salmonpy/distributed.py" "${ROOT}/salmonpy/material_response.py"
  "${ROOT}/salmonpy/DISTRIBUTED.md" "${ROOT}/salmonpy/DISTRIBUTED.en.md" DESTINATION salmonpy)
install(FILES "${ROOT}/tools/build_parallel_bridge.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/salmonpy/embedded_tddft.py" "${ROOT}/salmonpy/coupled.py"
  "${ROOT}/salmonpy/COUPLED.md" "${ROOT}/salmonpy/COUPLED.en.md" DESTINATION salmonpy)
install(FILES "${ROOT}/tools/build_embedded_tddft.py" DESTINATION salmonpy/_runtime_source/tools)
