set(STRESS_SOURCE_TESTS
  check_stress_occupation_weight.sh
  check_stress_ewald_g_sign.sh
  check_stress_fd_detail_levels.sh
  check_stress_virial_diagnostics.sh
  check_stress_decomp_tensor_table.sh
  check_stress_nl_factor_gate.sh
  check_stress_hartree_shadow.sh
  check_stress_local_sr_rs_shadow.sh
  check_stress_local_diag_gate.sh
  check_stress_local_grad_norm_gate.sh
  check_stress_local_grad_form_gate.sh
  check_stress_convention_a_sign_flip.sh
  check_stress_xc_split.sh
)

foreach(test_script IN LISTS STRESS_SOURCE_TESTS)
  get_filename_component(test_name "${test_script}" NAME_WE)
  add_test(
    NAME "${test_name}"
    COMMAND "${CMAKE_COMMAND}" -E env sh "${CMAKE_CURRENT_SOURCE_DIR}/${test_script}"
  )
endforeach()
