"""Freeze diagnostic densities before the terminal potential can change."""
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = (root / "src/gs/main_dft.f90").read_text().lower()
lcfo = (root / "src/gs/dc/lcfo.f90").read_text().lower()
begin = source.index("subroutine run_dg_hybrid_divided_ground_state_for_main")
end = source.index("end subroutine run_dg_hybrid_divided_ground_state_for_main", begin)
route = source[begin:end]
assert "retained_core_density" in lcfo, "conventional LCFO core-density output missing"
assert "subroutine build_retained_core_density" in lcfo
density = lcfo[lcfo.index("subroutine build_retained_core_density"):]
density = density[:density.index("end subroutine build_retained_core_density")]
assert "f_basis" in density and "coef_wf" in density and "retained_occupations" in density
assert "partition_weight" not in density and "spsi%rwf" not in density
assert "salmon_dg_density_diagnostic_prefix" in route
assert route.index("call update_dg_hybrid_divided_potential(initial_density") < route.index("call dc_lcfo(")
assert route.index("call dc_lcfo(") < route.index("call initialize_dg_hybrid_interface_continuation")
assert route.index("allocate(diagnostic_frozen_density,source=terminal_density_output)") < route.index(
    "call evaluate_dg_hybrid_terminal_total_energy")
assert route.index("call write_dg_density_diagnostic") > route.index("enddo terminal_lcfo_refinement")
print("frozen-potential density snapshot route: PASS")
