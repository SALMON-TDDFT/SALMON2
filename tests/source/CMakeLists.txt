set(STRESS_SOURCE_TESTS
  check_stress_occupation_weight.sh
  check_stress_ewald_g_sign.sh
  check_stress_ewald_rspace_sign.sh
  check_stress_output_control_keys.sh
  check_stress_virial_diagnostics.sh
  check_stress_decomp_tensor_table.sh
  check_stress_field_inventory_redesign.sh
  check_stress_output_level_profiles_redesign.sh
  check_stress_energy_schema_redesign.sh
  check_stress_rt_output_schema_redesign.sh
  check_stress_nl_factor_gate.sh
  check_stress_hartree_shadow.sh
  check_stress_local_sr_rs_shadow.sh
  check_stress_local_diag_gate.sh
  check_stress_local_grad_norm_gate.sh
  check_stress_local_grad_form_gate.sh
  check_stress_convention_a_sign_flip.sh
  check_stress_xc_split.sh
  check_stress_pz_nlcc_density_gate.sh
  check_prep_pp_dvg_ion_dg2_series.sh
  check_builtin_r2scan_registration.sh
  check_builtin_r2scan_loop_perf_structure.sh
  check_builtin_r2scan_operator_gate.sh
  check_builtin_r2scan_vtau_damping.sh
  check_builtin_r2scan_rt_ms_enablement.sh
  check_builtin_r2scan_nonorth_tau_operator.sh
  check_builtin_r2scan_nonorth_tau_operator_no_slice_temps.sh
  check_builtin_r2scan_primitive_smoke.sh
  check_bulk_si_rt_ms_stress_smoke_tool.sh
  check_nonorth_stress_keeps_orthogonal_branch.sh
  check_primitive_nonorth_stress_fd_tool.sh
  check_builtin_r2scan_stress_payload_consistency.sh
  check_builtin_r2scan_stress_detail_split.sh
  check_builtin_r2scan_stress_grho_payload_diag.sh
  check_builtin_r2scan_stress_gradient_trace_diag.sh
  check_builtin_r2scan_stress_vsigma_diag.sh
  check_builtin_r2scan_stress_guard.sh
  check_tau_mixing_input_keys.sh
  check_tau_mixing_broyden_pack.sh
  check_tau_mixing_xc_override.sh
)

foreach(test_script IN LISTS STRESS_SOURCE_TESTS)
  get_filename_component(test_name "${test_script}" NAME_WE)
  add_test(
    NAME "${test_name}"
    COMMAND "${CMAKE_COMMAND}" -E env sh "${CMAKE_CURRENT_SOURCE_DIR}/${test_script}"
  )
endforeach()
