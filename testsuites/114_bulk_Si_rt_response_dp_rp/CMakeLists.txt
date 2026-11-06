create_test()
# This response test reads the GS restart; order both calculations explicitly.
set_tests_properties(run_114_bulk_Si_rt_response_dp_rp PROPERTIES
  FIXTURES_REQUIRED "setup_114_bulk_Si_rt_response_dp_rp;bulk_si_dp_gs_ready"
  FIXTURES_SETUP run_114_bulk_Si_rt_response_dp_rp)
