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
add_library(salmonpy_ld SHARED "${ROOT}/native/lorentz_drude.f90" "${ROOT}/salmonpy/ld_c_api.f90" "${ROOT}/salmonpy/ld_handles.f90"
  "${ROOT}/native/sbe_types.f90" "${ROOT}/native/sbe_velocity_kernel.f90" "${ROOT}/salmonpy/sbe_c_api.f90" "${ROOT}/salmonpy/sbe_handles.f90")
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
# install(TARGETS) may change RPATH/install names. Hash the installed bytes,
# including staged DESTDIR installs; wheel repair recalculates again if needed.
install(CODE [[
  set(native_file "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/salmonpy/_native/$<TARGET_FILE_NAME:salmonpy_ld>")
  file(SHA256 "${native_file}" native_sha)
  file(WRITE "${native_file}.sha256" "${native_sha}\n")
]])
install(FILES "${ROOT}/salmonpy/__init__.py" "${ROOT}/salmonpy/ld.py"
  "${ROOT}/salmonpy/fields.py" "${ROOT}/salmonpy/geometry.py"
  "${ROOT}/salmonpy/salmon_vr.py" "${ROOT}/salmonpy/SALMON_VR.md" "${ROOT}/salmonpy/SALMON_VR.en.md" "${ROOT}/salmonpy/doctor.py" "${ROOT}/salmonpy/AGENTS.md" "${ROOT}/salmonpy/manual.py" "${ROOT}/salmonpy/manual.json" "${ROOT}/salmonpy/ase_calculator.py" "${ROOT}/salmonpy/dft.py" "${ROOT}/salmonpy/ase.py" "${ROOT}/salmonpy/bands.py" "${ROOT}/salmonpy/analysis.py" "${ROOT}/salmonpy/outputs.py" "${ROOT}/salmonpy/fdtd.py" "${ROOT}/salmonpy/crystal.py" "${ROOT}/salmonpy/sbe.py" "${ROOT}/salmonpy/tddft.py"
  "${ROOT}/salmonpy/ASE_CALCULATOR.md" "${ROOT}/salmonpy/ASE_CALCULATOR.en.md" "${ROOT}/salmonpy/DFT.md" "${ROOT}/salmonpy/ASE.md" "${ROOT}/salmonpy/ASE.en.md" "${ROOT}/salmonpy/BANDS.md" "${ROOT}/salmonpy/WORKFLOWS.md" "${ROOT}/salmonpy/OUTPUTS.md" "${ROOT}/salmonpy/PACKAGING.md" "${ROOT}/salmonpy/README.md" "${ROOT}/salmonpy/README.en.md" "${ROOT}/salmonpy/AI_USAGE.md" "${ROOT}/salmonpy/INPUT_CONTRACT.md"
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
  PATTERN "__pycache__" EXCLUDE PATTERN "*.pyc" EXCLUDE)
install(DIRECTORY "${ROOT}/native/" DESTINATION salmonpy/_runtime_source/native
  FILES_MATCHING PATTERN "*.f90")

install(FILES "${ROOT}/salmonpy/distributed.py" "${ROOT}/salmonpy/material_response.py"
  "${ROOT}/salmonpy/DISTRIBUTED.md" "${ROOT}/salmonpy/DISTRIBUTED.en.md" DESTINATION salmonpy)
install(FILES "${ROOT}/tools/build_parallel_bridge.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/salmonpy/maxwell.py" "${ROOT}/salmonpy/timing.py" "${ROOT}/salmonpy/FIELDS.md" "${ROOT}/salmonpy/FIELDS.en.md" DESTINATION salmonpy)
install(FILES "${ROOT}/salmonpy/embedded_tddft.py" "${ROOT}/salmonpy/coupled.py"
  "${ROOT}/salmonpy/COUPLED.md" "${ROOT}/salmonpy/COUPLED.en.md" DESTINATION salmonpy)
install(FILES "${ROOT}/tools/build_embedded_tddft.py" DESTINATION salmonpy/_runtime_source/tools)

# Fixed starter input, pseudopotential provenance and reference data.
install(DIRECTORY "${ROOT}/salmonpy/examples/si_intro_data/"
  DESTINATION salmonpy/examples/si_intro_data
  FILES_MATCHING PATTERN "*.json" PATTERN "*.dat" PATTERN "*.txt")
install(FILES "${ROOT}/salmonpy/SI_INTRO.md" "${ROOT}/salmonpy/SI_INTRO.en.md"
  DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/MEMORY.md" "${ROOT}/salmonpy/MEMORY.en.md" "${ROOT}/salmonpy/memory.py" "${ROOT}/salmonpy/_memory_advice.py" "${ROOT}/salmonpy/checkpointing.py" "${ROOT}/salmonpy/_native_checkpoint.py" "${ROOT}/salmonpy/CHECKPOINTING.md" "${ROOT}/salmonpy/CHECKPOINTING.en.md" DESTINATION salmonpy)

# Public Python integration contract and in-memory state snapshots.
install(FILES "${ROOT}/salmonpy/_version.py" "${ROOT}/salmonpy/errors.py"
  "${ROOT}/salmonpy/api.py" "${ROOT}/salmonpy/snapshots.py"
  "${ROOT}/salmonpy/API.md" "${ROOT}/salmonpy/API.en.md"
  "${ROOT}/salmonpy/SNAPSHOTS.md" "${ROOT}/salmonpy/SNAPSHOTS.en.md"
  DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/_tddft_checkpoint.py" "${ROOT}/salmonpy/_runtime_capabilities.py" "${ROOT}/salmonpy/_nlcc_reference.json"
  DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/materials.py" "${ROOT}/salmonpy/MATERIALS.md" "${ROOT}/salmonpy/MATERIALS.en.md" DESTINATION salmonpy)

# Preview public C material ABI; only declarations/examples, no implicit plugin build.
install(FILES "${ROOT}/salmonpy/c_material.py" "${ROOT}/salmonpy/C_ABI.md" "${ROOT}/salmonpy/C_ABI.en.md" DESTINATION salmonpy)
install(DIRECTORY "${ROOT}/salmonpy/include/" DESTINATION salmonpy/include
  FILES_MATCHING PATTERN "*.h" PATTERN "*.f90")
install(DIRECTORY "${ROOT}/salmonpy/examples/c_material/" DESTINATION salmonpy/examples/c_material
  FILES_MATCHING PATTERN "*.c" PATTERN "*.cpp" PATTERN "*.f90")
install(FILES "${ROOT}/tools/build_c_material_example.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/salmonpy/_sbe_handles.py" "${ROOT}/salmonpy/SBE.md" "${ROOT}/salmonpy/SBE.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/_ld_handles.py" "${ROOT}/salmonpy/LD.md" "${ROOT}/salmonpy/LD.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/gs_archive.py" "${ROOT}/salmonpy/GS_ARCHIVE.md" "${ROOT}/salmonpy/GS_ARCHIVE.en.md" DESTINATION salmonpy)
install(FILES "${ROOT}/salmonpy/include/salmonpy_gs_stream.f90" DESTINATION salmonpy/include)

install(FILES "${ROOT}/salmonpy/gs_native.py" DESTINATION salmonpy)
install(DIRECTORY "${ROOT}/salmonpy/native_adapters/" DESTINATION salmonpy/native_adapters FILES_MATCHING PATTERN "*.f90" PATTERN "*.h")
install(FILES "${ROOT}/tools/build_gs_restart_converter.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/salmonpy/_gs_export.py" DESTINATION salmonpy)

# Opt-in RT configuration extraction; default pinned runtime remains unchanged.
install(FILES "${ROOT}/tools/prepare_rt_step_config.py" "${ROOT}/tools/test_rt_step_config.py"
  "${ROOT}/tools/test_rt_step_config_runtime.py" DESTINATION salmonpy/_runtime_source/tools)
install(FILES "${ROOT}/tests/native/test_rt_step_config.f90" DESTINATION salmonpy/_runtime_source/tests/native)

# Separate opt-in fast-PR HSE profile; no full solver is bundled.
install(FILES "${ROOT}/tools/prepare_hse_rt_state.py" "${ROOT}/tools/test_hse_rt_state.py"
  "${ROOT}/tools/test_hse_rt_state_runtime.py" DESTINATION salmonpy/_runtime_source/tools)
install(FILES "${ROOT}/tests/native/test_hse_rt_state.f90" DESTINATION salmonpy/_runtime_source/tests/native)

install(FILES "${ROOT}/tools/test_pp_grid_config.py" "${ROOT}/tools/test_kpoint_config.py" "${ROOT}/tools/test_pp_input_config.py" "${ROOT}/tools/test_xc_input_config.py" "${ROOT}/tools/test_pre_scf_state.py" "${ROOT}/tools/test_pre_scf_runtime.py" "${ROOT}/tools/test_gs_lifecycle.py" "${ROOT}/tools/test_gs_occupation_runtime.py" "${ROOT}/tools/test_occupation_inputs.py" "${ROOT}/tools/test_cell_config.py" "${ROOT}/tools/test_grid_config.py" "${ROOT}/tools/test_owned_workspaces.py" "${ROOT}/tools/test_tddft_contexts.py" "${ROOT}/tools/test_hse_explicit_state.py" DESTINATION salmonpy/_runtime_source/tools)
install(FILES "${ROOT}/tests/native/test_hse_explicit_state.f90" DESTINATION salmonpy/_runtime_source/tests/native)

install(FILES "${ROOT}/tests/native/test_pp_grid_config.f90" "${ROOT}/tests/native/test_kpoint_config.f90" "${ROOT}/tests/native/test_pp_input_config.f90" "${ROOT}/tests/native/test_xc_input_config.f90" "${ROOT}/tests/native/test_pre_scf_state.f90" "${ROOT}/tests/native/test_gs_lifecycle.f90" "${ROOT}/tests/native/test_occupation_inputs.f90" "${ROOT}/tests/native/test_dip_occupation_inputs.f90" "${ROOT}/tests/native/test_cell_config.f90" "${ROOT}/tests/native/test_grid_config.f90" "${ROOT}/tests/native/test_owned_workspaces.f90" "${ROOT}/tests/native/test_tddft_contexts.f90" "${ROOT}/tests/native/test_hse_native_resources.f90" "${ROOT}/tests/native/test_hse_resource_context.f90" DESTINATION salmonpy/_runtime_source/tests/native)

install(FILES "${ROOT}/tools/runtime_sources.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/salmonpy/config.py" "${ROOT}/salmonpy/CONFIG.md" "${ROOT}/salmonpy/CONFIG.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/preparation.py" "${ROOT}/salmonpy/PREPARATION.md" "${ROOT}/salmonpy/PREPARATION.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/structure.py" "${ROOT}/salmonpy/STRUCTURE.md" "${ROOT}/salmonpy/STRUCTURE.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/builtin_materials.py" "${ROOT}/salmonpy/BUILTIN_MATERIALS.md" "${ROOT}/salmonpy/BUILTIN_MATERIALS.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/DFT_CONFIG.md" "${ROOT}/salmonpy/DFT_CONFIG.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/initial_states.py" "${ROOT}/salmonpy/INITIAL_STATES.md" "${ROOT}/salmonpy/INITIAL_STATES.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/point_runs.py" "${ROOT}/salmonpy/POINT_RUNS.md" "${ROOT}/salmonpy/POINT_RUNS.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/MPI_IDENTITY.md" "${ROOT}/salmonpy/MPI_IDENTITY.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/RECORDED_MATERIALS.md" "${ROOT}/salmonpy/RECORDED_MATERIALS.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/CONFIGURED_NATIVE.md" "${ROOT}/salmonpy/CONFIGURED_NATIVE.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/_point_checkpoint.py" "${ROOT}/salmonpy/POINT_CHECKPOINTS.md" "${ROOT}/salmonpy/POINT_CHECKPOINTS.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/field_input.py" "${ROOT}/salmonpy/FIELD_INPUT.md" "${ROOT}/salmonpy/FIELD_INPUT.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/dipole_source.py" "${ROOT}/salmonpy/DIPOLE_RADIATION.md" "${ROOT}/salmonpy/DIPOLE_RADIATION.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/isolated_tddft.py" "${ROOT}/salmonpy/ISOLATED_TDDFT.md" "${ROOT}/salmonpy/ISOLATED_TDDFT.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/field_sampling.py" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/layers.py" DESTINATION salmonpy)
install(FILES "${ROOT}/salmonpy/examples/mixed_layers.json" DESTINATION salmonpy/examples)

# Fast-PR HSE direct numerical checkpoint regression, explicitly built only.
install(FILES "${ROOT}/tools/test_hse_direct_checkpoint.py" "${ROOT}/tools/test_fast_hse_mpi_checkpoint.py" DESTINATION salmonpy/_runtime_source/tools)
install(FILES "${ROOT}/tests/native/test_hse_direct_checkpoint.f90" DESTINATION salmonpy/_runtime_source/tests/native)

# Explicit process-isolated fast-PR HSE material adapter and worker sources.
install(FILES "${ROOT}/salmonpy/fast_hse.py" "${ROOT}/salmonpy/fast_hse_mpi.py" DESTINATION salmonpy)
install(FILES "${ROOT}/salmonpy/examples/mixed_hse_layers.json" DESTINATION salmonpy/examples)
install(FILES "${ROOT}/tools/build_fast_hse_worker.py" "${ROOT}/tools/build_fast_hse_mpi.py" "${ROOT}/tools/prepare_fast_hse_coupled_fixture.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/tools/test_native_tddft_admission_mpi.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/tools/test_mixed_scaling.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/tools/test_tddft_deep_recovery_mpi.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/tools/test_fast_hse_deep_recovery_mpi.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/tools/test_eh_deep_recovery_mpi.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/salmonpy/stepping.py" "${ROOT}/salmonpy/STEPPING.md" "${ROOT}/salmonpy/STEPPING.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/native_steps.py" "${ROOT}/salmonpy/NATIVE_STEPPING.md" "${ROOT}/salmonpy/NATIVE_STEPPING.en.md" DESTINATION salmonpy)
install(FILES "${ROOT}/salmonpy/examples/mixed_tddft_hse_layers.json" DESTINATION salmonpy/examples)
install(FILES "${ROOT}/tools/test_native_protected_step_mpi.py" "${ROOT}/tools/test_native_step_admission_mpi.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/salmonpy/examples/mixed_multipoint_layers.json" DESTINATION salmonpy/examples)

install(FILES "${ROOT}/tests/native/test_ion_initialization.f90" DESTINATION salmonpy/_runtime_source/tests/native)

install(FILES "${ROOT}/tests/native/test_orbital_initialization.f90" DESTINATION salmonpy/_runtime_source/tests/native)

install(FILES "${ROOT}/tests/native/test_hse_symmetry_resources.f90" DESTINATION salmonpy/_runtime_source/tests/native)
install(FILES "${ROOT}/tests/native/test_hse_symmetry_field.f90" DESTINATION salmonpy/_runtime_source/tests/native)
install(FILES "${ROOT}/tests/native/test_hse_symmetry_checkpoint.f90" DESTINATION salmonpy/_runtime_source/tests/native)
install(FILES "${ROOT}/tests/native/test_hse_symmetry_tddft.f90" DESTINATION salmonpy/_runtime_source/tests/native)
install(FILES "${ROOT}/tools/prepare_fast_hse_symmetry_fixture.py" "${ROOT}/tools/test_fast_hse_symmetry_tddft.py" DESTINATION salmonpy/_runtime_source/tools)
install(FILES "${ROOT}/tools/test_fast_hse_symmetry_mpi.py" "${ROOT}/tools/probe_hse_low_gradient.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/salmonpy/FAST_HSE_SYMMETRY.md" "${ROOT}/salmonpy/FAST_HSE_SYMMETRY.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/tools/test_hse_low_gradient.py" "${ROOT}/tools/test_fast_hse_low_gradient_mpi.py" DESTINATION salmonpy/_runtime_source/tools)
install(FILES "${ROOT}/tests/native/test_hse_low_gradient.f90" DESTINATION salmonpy/_runtime_source/tests/native)

install(FILES "${ROOT}/salmonpy/_hse_geometry.py" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/_hse_assets.py" DESTINATION salmonpy)

install(FILES "${ROOT}/tools/test_fast_hse_general_cell_mpi.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/tools/test_fast_hse_symmetry_history_mpi.py" DESTINATION salmonpy/_runtime_source/tools)
install(FILES "${ROOT}/tests/native/test_hse_context_owners.f90" "${ROOT}/tests/native/test_step_output_owners.f90" DESTINATION salmonpy/_runtime_source/tests/native)
install(FILES "${ROOT}/tools/test_fast_hse_context_owners_mpi.py" DESTINATION salmonpy/_runtime_source/tools)
install(FILES "${ROOT}/tests/native/test_tddft_output_close.f90" DESTINATION salmonpy/_runtime_source/tests/native)

install(FILES "${ROOT}/tools/test_fast_hse_handles_mpi.py" DESTINATION salmonpy/_runtime_source/tools)
install(FILES "${ROOT}/salmonpy/examples/native_fast_hse_handles.c" DESTINATION salmonpy/examples)

install(FILES "${ROOT}/salmonpy/NATIVE_HANDLES.md" "${ROOT}/salmonpy/NATIVE_HANDLES.en.md" DESTINATION salmonpy)

install(FILES "${ROOT}/salmonpy/_shared_hse.py" DESTINATION salmonpy)
install(FILES "${ROOT}/tools/test_shared_hse_mpi.py" DESTINATION salmonpy/_runtime_source/tools)

install(FILES "${ROOT}/salmonpy/SHARED_HSE.md" "${ROOT}/salmonpy/SHARED_HSE.en.md" DESTINATION salmonpy)
