#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src/gs/dc/dg_wpw_matrix_free_scf.f90"
MAIN = ROOT / "src/gs/main_dft.f90"
LCFO = ROOT / "src/gs/dc/lcfo_flux.f90"
GLOBAL = ROOT / "src/io/salmon_global.f90"
INPUT = ROOT / "src/io/inputoutput.f90"

assert SOURCE.exists(), "missing bounded matrix-free DG-DC consumer"
source = SOURCE.read_text().lower()
main = MAIN.read_text().lower()
lcfo = LCFO.read_text().lower()
global_text = GLOBAL.read_text().lower()
input_text = INPUT.read_text().lower()
controls = global_text + input_text

precondition_control = "yn_dg_wpw_preconditioner"
assert precondition_control in global_text, "missing explicit WPW preconditioner comparison control"
assert precondition_control in input_text, "WPW preconditioner control is absent from input plumbing"
assert f"{precondition_control} = 'y'" in input_text, "current preconditioned route must remain the default"
for required in ("call comm_bcast(yn_dg_wpw_preconditioner", "call yn_argument_check(yn_dg_wpw_preconditioner)"):
    assert required in input_text, f"missing WPW preconditioner input contract: {required}"
assert lcfo.count("if(yn_dg_wpw_preconditioner=='y')then") >= 2, \
    "fixed-H and continuation routes must independently select the optional preconditioner"
for routine in ("run_wpw_fixed_h_relaxation", "continue_wpw_fixed_h_interface"):
    start = lcfo.index(f"subroutine {routine}")
    end = lcfo.index("end subroutine", start)
    body = lcfo[start:end]
    branch = body.index("if(yn_dg_wpw_preconditioner=='y')then")
    fallback = body.index("else", branch)
    closing = body.index("endif", fallback)
    assert "wpw_precondition" in body[branch:fallback]
    assert "wpw_precondition" not in body[fallback:closing]

for token in (
    "run_dg_wpw_matrix_free_scf",
    "apply_h_batch",
    "apply_s_batch",
    "global_gram_batch",
    "n_occ + extra_states",
    "density_residual",
    "potential_residual",
    "energy_residual",
    "projector_residual",
    "generalized_residual",
    "metric_orthonormality",
    "charge_error",
    "gap",
    "fixed_point",
):
    assert token in source, f"missing matrix-free SCF contract: {token}"

for forbidden in (
    "allocate(h_global",
    "allocate(s_global",
    "allocate(density_matrix",
    "run_dg_wpw_fixed_scf",
):
    assert forbidden not in source + main, f"dense production path detected: {forbidden}"

for stage in (
    "initial_apply_begin", "initial_apply_end", "initial_eigh_end",
    "expanded_apply_end", "reduced_eigh_end", "inner_end",
):
    assert f"call trace_solve_window('{stage}'" in source, \
        f"missing solve-window boundary diagnostic: {stage}"
assert "dg_wpw_preconditioner" in source and "present(precondition)" in source, \
    "matrix-free solver must accept an optional bounded preconditioner callback"
solve_start = source.index("subroutine solve_window")
solve_end = source.index("end subroutine", solve_start)
solve_body = source[solve_start:solve_end]
assert "nw,np,nocc,nretain" in solve_body[:500], "solve_window must know the occupied/extra partition"
assert "window_state_residual_iteration(inner)" in solve_body, "state diagnostics need bounded selected iterations"
assert "call gram(preconditioned,preconditioned" in solve_body, \
    "preconditioner response must use the production global Gram"
assert "occupied_preconditioned=" in solve_body and "extra_preconditioned=" in solve_body, \
    "diagnostic must print preconditioned norms as well as ratios"
state_diag_pos = solve_body.index("[dg-wpw-window-state-residual]")
convergence_pos = solve_body.index("if(residual<tol.and.orth<tol)return")
assert convergence_pos < state_diag_pos, "state diagnostics must not replace the existing convergence decision"
assert "residual=max(occupied_max,extra_max)" not in solve_body, \
    "split diagnostics must not feed back into convergence"
metric_mode_pos = solve_body.index("[dg-wpw-search-metric-mode]")
assert "call gram(z,r" in solve_body[:metric_mode_pos], "search metric diagnostic must assemble Z^H R"
assert "call dg_metric_mode_residual_split" in solve_body[:metric_mode_pos]
assert metric_mode_pos < solve_body.index("call dg_reduced_generalized_eigh"), \
    "metric-mode diagnostic must observe, not replace, the production reduced solve"

for token in (
    "yn_dg_wpw_production",
    "dg_wpw_extra_states",
    "dg_wpw_gap_threshold",
    "dg_wpw_metric_cutoff",
    "dg_wpw_scf_mix",
    "dg_wpw_scf_max_iter",
    "dg_wpw_scf_residual_tolerance",
):
    assert token in controls, f"missing DG-DC namelist control: {token}"

mod = Path("/tmp/test_dg_wpw_matrix_free_scf_mod")
mod.mkdir(exist_ok=True)
exe = Path("/tmp/test_dg_wpw_matrix_free_scf_mpi")
fc = os.environ.get("MPIFC", "mpifort")
subprocess.run([fc,"-J",str(mod),"-I",str(mod),
    str(ROOT/"src/common/dg_generalized_algebra.f90"),
    str(ROOT/"src/common/dg_wpw_matrix_free_operator.f90"),str(SOURCE),
    str(ROOT/"tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90"),"-llapack","-lblas","-o",str(exe)],check=True,cwd="/tmp")
subprocess.run(["mpirun","-np","2",str(exe)],check=True,cwd=ROOT)
assert "call dc_lcfo_flux" in main, "main_dft does not reach the LCFO production owner"
assert "call run_wpw_production_scf" in lcfo, "lcfo_flux does not reach the stepwise production consumer"
print("PASS bounded matrix-free DG-DC production route contract")
