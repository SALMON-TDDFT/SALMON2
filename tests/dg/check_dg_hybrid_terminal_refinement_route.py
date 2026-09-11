#!/usr/bin/env python3
"""Route boundary for local Hybrid continuation and terminal LCFO refinement."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()
INPUT = (ROOT / "src/io/inputoutput.f90").read_text(errors="replace").lower()


def subroutine(name: str) -> str:
    begin = f"subroutine {name}"
    end = f"end subroutine {name}"
    assert begin in SOURCE and end in SOURCE, f"missing {name}"
    return SOURCE[SOURCE.index(begin) : SOURCE.index(end, SOURCE.index(begin))]


alias = subroutine("run_dg_hybrid_continuation_ground_state_for_main")
assert alias.count("call run_dg_hybrid_divided_ground_state_for_main") == 1, (
    "the continuation input must alias the one-fragment-per-rank divided route"
)
assert "run_dg_overlapping_wannier_ground_state_for_main" not in alias, (
    "the continuation input still enters the legacy whole-system stage loop"
)

dispatch_begin = SOURCE.index("if(yn_dg_hybrid_divided_scf == 'y') then")
dispatch_end = SOURCE.index("end if", dispatch_begin)
dispatch = SOURCE[dispatch_begin:dispatch_end]
assert "call run_dg_hybrid_divided_ground_state_for_main" in dispatch
assert "call run_dg_hybrid_continuation_ground_state_for_main" in dispatch
assert dispatch.index("yn_dg_hybrid_divided_scf") < dispatch.index(
    "yn_dg_hybrid_continuation_scf"
), "the explicit divided selector must retain priority"

divided = subroutine("run_dg_hybrid_divided_ground_state_for_main")
local_begin = divided.index("call initialize_dg_hybrid_interface_continuation")
terminal_begin = divided.index("call validate_dg_hybrid_schwarz_dynamic_receipt")
assert local_begin < terminal_begin
local_phase = divided[local_begin:terminal_begin]
terminal_phase = divided[terminal_begin:]

assert "call solve_dg_hybrid_schwarz_fragments" in local_phase
assert "solve_dg_hybrid_generalized" not in local_phase, (
    "fragment-local continuation must not perform a complete eigensolve"
)
assert "rho=" not in local_phase and "initial_density=" not in local_phase, (
    "the local continuation must not replace the immutable DC seed density"
)
assert "interface_continuation%lambda" in local_phase
assert "bounded_interface_scale" in local_phase

assert terminal_phase.count("call solve_dg_hybrid_generalized_once_and_publish") == 1
assert "final_hrows=bounded_fixed_payload%kinetic_rows+bounded_fixed_payload%nonlocal_rows+&" in terminal_phase
assert "bounded_fixed_payload%interface_rows+final_local_potential_rows" in terminal_phase
assert "final_srows=bounded_fixed_payload%metric_rows" in terminal_phase
assert "call initialize_dg_hybrid_terminal_refinement" in terminal_phase, (
    "terminal LCFO refinement policy is not initialized"
)
assert "full_from_start=.true." in local_phase, (
    "the fragment-local phase still ramps the DG interface instead of applying it fully from step one"
)
assert "call initialize_dg_hybrid_terminal_operator_guard" in terminal_phase, (
    "terminal LCFO immutable components are not fingerprinted"
)
assert "terminal_lcfo_refinement: do" in terminal_phase
refinement_loop = terminal_phase[terminal_phase.index("terminal_lcfo_refinement: do") :]
assert refinement_loop.count("call solve_dg_hybrid_generalized_once_and_publish") == 1, (
    "one syntactic eigensolve call must serve all one-to-four terminal solves"
)
assert "call reconstruct_dg_hybrid_terminal_density" in refinement_loop
assert "call mix_dg_overlapping_wannier_density_history" in refinement_loop
assert "call observe_dg_hybrid_terminal_refinement" in refinement_loop
assert refinement_loop.count("call validate_dg_hybrid_terminal_operator_guard") >= 2, (
    "immutable components must be checked before and after each terminal solve"
)
assert "if(.not.terminal_request_another)exit terminal_lcfo_refinement" in refinement_loop
assert "dg_dc_gs_final_density_tolerance" in terminal_phase
assert "ow_hybrid_divided_threshold" in terminal_phase, (
    "terminal energy convergence must reuse the conventional DC threshold"
)
assert SOURCE.count("[dg-hybrid-refinement-warning] maximum additional lcfo solves exhausted; publishing last finite valid state") == 1, (
    "terminal exhaustion warning must be emitted exactly once"
)
for field in (
    "final_eigensolve_count",
    "additional_refinement_count",
    "refinement_converged",
    "refinement_exhausted",
    "terminal_density_change",
    "terminal_energy_change",
):
    assert f"ow_hybrid_ground_state%{field}" in terminal_phase, (
        f"terminal ground-state audit omits {field}"
    )

fragment_solver = subroutine("solve_dg_hybrid_schwarz_fragments")
assert "solve_dg_hybrid_generalized" not in fragment_solver
assert "dg_hybrid_fragment_cg_steps" in fragment_solver
assert "local_iterations<=dg_hybrid_fragment_cg_steps" in fragment_solver, (
    "fragment-local CG work must obey the configured small-step cap"
)
assert "dg_hybrid_fragment_cg_steps = 3" in INPUT, (
    "fragment-local CG must retain the three-step production default"
)

print("terminal LCFO/local continuation route boundary: PASS")
