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

complement_control = "yn_dg_wpw_s_orthogonal_pw"
assert complement_control in global_text, "missing S-orthogonal PW comparison control"
assert complement_control in input_text, "S-orthogonal PW control is absent from input plumbing"
assert f"{complement_control} = 'n'" in input_text, "S-orthogonal PW comparison must remain default-off"
for required in (
    "call comm_bcast(yn_dg_wpw_s_orthogonal_pw",
    "'yn_dg_wpw_s_orthogonal_pw', yn_dg_wpw_s_orthogonal_pw",
    "call yn_argument_check(yn_dg_wpw_s_orthogonal_pw)",
):
    assert required in input_text, f"missing S-orthogonal PW input contract: {required}"
assert "yn_dg_wpw_s_orthogonal_pw=='y'.and.yn_dg_wpw_fixed_h_relaxation/='y'" in \
       input_text.replace(" ", ""), "S-orthogonal PW comparison must be fixed-H-only"
assert "initialize_dg_wpw_s_orthogonal_complement" in lcfo
assert "wpw_apply_h_complement" in lcfo and "wpw_apply_s_complement" in lcfo
assert "wpw_fixed_apply_h=>wpw_apply_h_complement" in lcfo.replace(" ", "")
assert "wpw_fixed_apply_h=>wpw_apply_h" in lcfo.replace(" ", "")
assert "[dg-wpw-pw-complement] mode=s-orthogonal fingerprint=" in lcfo
for field in ("solve_residual=", "cross_metric_defect=", "numerical_p_rank="):
    assert field in lcfo, f"missing complement log field: {field}"
assert "map_dg_wpw_original_to_complement" in lcfo
assert lcfo.count("map_dg_wpw_complement_to_original") >= 2, \
    "density and checkpoint boundaries must restore original coordinates"

precondition_control = "yn_dg_wpw_preconditioner"
metric_precondition_control = "yn_dg_wpw_metric_preconditioner"
assert precondition_control in global_text, "missing explicit WPW preconditioner comparison control"
assert precondition_control in input_text, "WPW preconditioner control is absent from input plumbing"
assert f"{precondition_control} = 'y'" in input_text, "current preconditioned route must remain the default"
for required in ("call comm_bcast(yn_dg_wpw_preconditioner", "call yn_argument_check(yn_dg_wpw_preconditioner)"):
    assert required in input_text, f"missing WPW preconditioner input contract: {required}"
assert metric_precondition_control in global_text, "missing metric-block correction comparison control"
assert metric_precondition_control in input_text, "metric-block correction control is absent from input plumbing"
assert f"{metric_precondition_control} = 'n'" in input_text, \
    "metric-block correction must remain default-off"
for required in (
    "call comm_bcast(yn_dg_wpw_metric_preconditioner",
    "'yn_dg_wpw_metric_preconditioner', yn_dg_wpw_metric_preconditioner",
    "call yn_argument_check(yn_dg_wpw_metric_preconditioner)",
):
    assert required in input_text, f"missing metric-block correction input contract: {required}"
assert "yn_dg_wpw_preconditioner=='y'.and.yn_dg_wpw_metric_preconditioner=='y'" in \
       input_text.replace(" ", ""), "diagonal and metric-block corrections must be mutually exclusive"
for routine in ("run_wpw_fixed_h_relaxation", "continue_wpw_fixed_h_interface"):
    start = lcfo.index(f"subroutine {routine}")
    end = lcfo.index("end subroutine", start)
    body = lcfo[start:end]
    metric_branch = body.index("if(yn_dg_wpw_metric_preconditioner=='y')then")
    diagonal_branch = body.index("elseif(yn_dg_wpw_preconditioner=='y')then", metric_branch)
    fallback = body.index("\n          else\n", diagonal_branch)
    closing = body.index("endif", fallback)
    assert "wpw_metric_precondition" in body[metric_branch:diagonal_branch]
    assert "wpw_precondition" in body[diagonal_branch:fallback]
    assert "precondition=" not in body[fallback:closing]
    assert body[metric_branch:closing].count("retain_search_history=yn_dg_wpw_search_history=='y'") == 3

algebra_start = lcfo.index("subroutine wpw_algebra_step")
algebra_end = lcfo.index("end subroutine wpw_algebra_step", algebra_start)
algebra_body = lcfo[algebra_start:algebra_end]
assert "precondition=" not in algebra_body, "normal production algebra must remain callback-free"
assert "retain_search_history=yn_dg_wpw_search_history=='y'" in algebra_body
assert "[dg-wpw-correction-mode]" in lcfo, "effective fixed-H correction mode must be logged"

history_control = "yn_dg_wpw_search_history"
assert history_control in global_text, "missing explicit WPW search-history comparison control"
assert history_control in input_text, "WPW search-history control is absent from input plumbing"
assert f"{history_control} = 'y'" in input_text, "current search-history route must remain the default"
for required in (
    "call comm_bcast(yn_dg_wpw_search_history",
    "'yn_dg_wpw_search_history', yn_dg_wpw_search_history",
    "call yn_argument_check(yn_dg_wpw_search_history)",
):
    assert required in input_text, f"missing WPW search-history input contract: {required}"
assert lcfo.count("retain_search_history=yn_dg_wpw_search_history=='y'") >= 3, \
    "all production algebra routes must propagate the search-history choice"

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
assert "retain_search_history" in solve_body[:600], "solve_window must accept the restart comparison flag"
assert "keep_search=.true." in solve_body and "present(retain_search_history)" in solve_body, \
    "omitting the comparison flag must retain the current search-history behavior"
search_update = solve_body.index("search=matmul(preconditioned")
search_branch = solve_body.rfind("if(keep_search)then", 0, search_update)
assert search_branch >= 0 and "else;search=(0d0,0d0);endif" in solve_body[search_update:], \
    "explicit restart mode must skip the discarded history update and clear only the search block"
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
    "ritz_consistency_comparison_iteration(inner)",
    "ritz_consistency_arm_iteration(inner)",
    "pending_ritz_residual",
    "direct_ritz_residual-pending_ritz_residual",
    "matmul(reduced_h,reduced_c(:,i))",
    "matmul(reduced_s,reduced_c(:,i))",
    "call mpi_allreduce(local_bad,global_bad",
    "call gram(q,post_ritz_sq",
    "[dg-wpw-ritz-consistency]",
):
    assert token in solve_body, f"missing Ritz-boundary diagnostic contract: {token}"
assert "max(direct_norm,post_norm)" in source and "direct_norm==0d0.and.post_norm==0d0" in source, \
    "Ritz relative vector defects need symmetric zero-safe normalization"
assert "inner==31.or.inner==95.or.inner==159" in source
assert "inner==32.or.inner==96.or.inner==160" in source
ritz_pos = solve_body.index("[dg-wpw-ritz-consistency]")
assert convergence_pos < solve_body.index("if(ritz_consistency_arm_iteration"), \
    "Ritz pending state must only be armed after the unchanged convergence test fails"
assert solve_body.index("call trace_solve_window('initial_apply_end'") < ritz_pos, \
    "Ritz comparison must reuse the next production direct H/S application"
compare_start = solve_body.index("if(pending_ritz_comparison.and.ritz_consistency_comparison_iteration(inner))then")
compare_end = solve_body.index("call gram(q,hq", compare_start)
compare_body = solve_body[compare_start:compare_end]
assert compare_body.count("call mpi_allreduce(local_bad,global_bad") >= 2
assert compare_body.count("pending_ritz_comparison=.false.") >= 2
arm_start = solve_body.index("if(ritz_consistency_arm_iteration(inner))then")
arm_end = solve_body.index("call trace_solve_window('inner_end'", arm_start)
arm_body = solve_body[arm_start:arm_end]
assert arm_body.count("call mpi_allreduce(local_bad,global_bad") >= 2
assert "pending_ritz_comparison=.false." in arm_body and "pending_ritz_comparison=.true." in arm_body
assert "pending_inner+1==inner" in compare_body

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
run = subprocess.run(["mpirun","-np","2",str(exe)],check=True,cwd=ROOT,
    stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
ritz_lines = [line for line in run.stdout.splitlines() if "[DG-WPW-RITZ-CONSISTENCY]" in line]
print("\n".join(ritz_lines))
for post_inner,direct_inner in ((31,32),(95,96),(159,160)):
    line = None
    values = None
    for candidate in ritz_lines:
        fields = candidate.replace("=", " ").split()
        parsed = {fields[i]: float(fields[i+1]) for i in range(len(fields)-1)
            if fields[i] in {"post_inner", "direct_inner", "relative_occupied", "relative_extra", "post_metric_orth"}}
        if int(parsed["post_inner"]) == post_inner and int(parsed["direct_inner"]) == direct_inner:
            line = candidate; values = parsed; break
    assert line is not None, f"missing runtime Ritz comparison {post_inner}->{direct_inner}"
    assert values["relative_occupied"] < 1e-8 and values["relative_extra"] < 1e-8, line
    assert values["post_metric_orth"] < 1e-8, line
assert "call dc_lcfo_flux" in main, "main_dft does not reach the LCFO production owner"
assert "call run_wpw_production_scf" in lcfo, "lcfo_flux does not reach the stepwise production consumer"
print("PASS bounded matrix-free DG-DC production route contract")
