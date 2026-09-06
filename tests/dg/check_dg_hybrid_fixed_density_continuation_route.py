#!/usr/bin/env python3
"""Source contract for the fixed-DC-density DG interface continuation route."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/main_dft.f90").read_text().lower()
source = re.sub(r"!.*", "", source).replace("&", "")

entry_name = "run_dg_hybrid_divided_ground_state_for_main"
match = re.search(
    rf"\bsubroutine\s+{entry_name}\b(.*?)\bend subroutine\s+{entry_name}",
    source,
    re.S,
)
assert match, "missing divided Hybrid production entry"
route = match.group(1)
compact = re.sub(r"\s+", "", route)

assert "nproc==dc%n_frag" in compact and "dc%isize_frag==1" in compact, (
    "fixed-density continuation must retain the exact rank-fragment guard"
)
assert "prepare_dg_hybrid_divided_dc_controls" in route, (
    "fixed-density continuation must restore the compatible ordinary-DC density"
)
assert re.search(
    r"call\s+initialize_dg_hybrid_interface_continuation\s*\(.*?mixing%mixrate",
    route,
    re.S,
), "production route must initialize interface continuation from mixing%mixrate"
assert "dg_hybrid_max_interface_points" in route, (
    "production continuation needs an explicit finite point budget"
)
assert re.search(
    r"interface_continuation%rate\s*>=\s*1d0\s*/\s*real\s*\(\s*"
    r"dg_hybrid_max_interface_points\s*-\s*2_int64",
    route,
), "rate must be range-checked before evaluating its reciprocal"
rate_guard = compact.index("interface_continuation%rate>=1d0/real")
rate_reciprocal = compact.index("1d0/interface_continuation%rate")
assert rate_guard < rate_reciprocal, "unsafe continuation-rate reciprocal precedes its guard"
assert "ceiling(1d0/interface_continuation%rate,kind=int64)+2_int64" in compact, (
    "continuation point count must use an explicitly wide integer kind"
)

potential_calls = re.findall(r"call\s+update_dg_hybrid_divided_potential\s*\(", route)
assert len(potential_calls) == 1, "restored DC density must update the potential exactly once"
assert re.search(r"call\s+update_dg_hybrid_divided_potential\s*\(\s*initial_density", route), (
    "the one potential update must use the immutable restored DC density"
)
assert route.index("call update_dg_hybrid_divided_potential") < route.index(
    "call initialize_dg_hybrid_interface_continuation"
), "fixed DC potential must be established before the continuation loop"

for obsolete_import in (
    "dg_hybrid_divided_scf",
    "dg_hybrid_divided_mixing",
):
    assert not re.search(r"\buse\s+" + obsolete_import + r"\b", source), (
        f"obsolete density-mixed production import remains: {obsolete_import}"
    )
for obsolete_state in (
    "divided_mixing_state",
    "divided_mixing_basis_generation",
    "divided_mixing_inventory_fingerprint",
    "divided_mixing_rollback_pending",
):
    assert obsolete_state not in source, (
        f"obsolete density-mixing production state remains: {obsolete_state}"
    )
for obsolete_callback in (
    "assemble_dg_hybrid_schwarz_core_density",
    "gather_dg_hybrid_divided_core_density",
    "mix_dg_hybrid_divided_density",
):
    assert not re.search(r"\bsubroutine\s+" + obsolete_callback + r"\b", source), (
        f"obsolete density-feedback callback remains: {obsolete_callback}"
    )

for forbidden in (
    "run_dg_hybrid_divided_scf",
    "assemble_dg_hybrid_schwarz_core_density",
    "mix_dg_hybrid_divided_density",
    "prepare_dg_hybrid_divided_mixing",
    "accept_dg_hybrid_divided_mixing",
    "pulay",
    "broyden",
):
    assert not re.search(r"\bcall\s+" + forbidden + r"\b", route), (
        f"density-feedback operation remains in fixed-density route: {forbidden}"
    )
for lifecycle in (
    "divided_mixing_basis_generation=",
    "divided_mixing_inventory_fingerprint=",
    "divided_mixing_rollback_pending=",
):
    assert lifecycle not in compact, (
        f"fixed-density production entry mutates density-mixing lifecycle: {lifecycle}"
    )

assert "dowhile(.not.interface_continuation%finished)" in compact, (
    "production route must solve every lambda point through terminal acceptance"
)
assert "bounded_interface_scale=interface_continuation%lambda" in compact, (
    "callbacks must receive the explicit live continuation lambda"
)
assert re.search(r"call\s+solve_dg_hybrid_schwarz_fragments\s*\(", route), (
    "each lambda point must run one capped Schwarz epoch"
)
assert "dg_hybrid_fragment_cg_steps" in source.split(
    "subroutine solve_dg_hybrid_schwarz_fragments", 1
)[1].split("end subroutine solve_dg_hybrid_schwarz_fragments", 1)[0], (
    "fixed-density Schwarz point must honor the configured CG cap"
)
assert "300d0" in source.split(
    "subroutine solve_dg_hybrid_schwarz_fragments", 1
)[1].split("end subroutine solve_dg_hybrid_schwarz_fragments", 1)[0], (
    "fixed-density continuation must assign common 300 K occupations"
)
solver = source.split("subroutine solve_dg_hybrid_schwarz_fragments", 1)[1].split(
    "end subroutine solve_dg_hybrid_schwarz_fragments", 1
)[0]
solver_compact = re.sub(r"\s+", "", solver)
assert "callback_ok=callback_ok.and..not.rolled_back" in solver_compact, (
    "a Schwarz rollback must be promoted to collective point failure"
)
collective_boundaries = (
    "extract_dg_hybrid_core_local_potential",
    "assemble_dg_hybrid_local_potential_rows",
    "ow_fingerprint_distributed_matrix",
    "assemble_dg_hybrid_schwarz_local_preconditioner_blocks",
    "advance_dg_hybrid_schwarz_epoch",
    "assign_dg_hybrid_schwarz_occupations",
)
for before, after in zip(collective_boundaries, collective_boundaries[1:]):
    segment = solver[solver.index(before) : solver.index(after)]
    assert "comm_logical_and" in segment, (
        f"missing collective failure gate between {before} and {after}"
    )
assert "allocate(core_potential(size(bounded_core_ids)),stat=status_local)" in solver_compact
assert solver_compact.index("stat=status_local") < solver_compact.index(
    "extract_dg_hybrid_core_local_potential"
)
for lifecycle in (
    "divided_mixing_inventory_fingerprint=",
    "divided_mixing_rollback_pending=",
):
    assert lifecycle not in re.sub(r"\s+", "", solver), (
        f"live Schwarz callback mutates density-mixing lifecycle: {lifecycle}"
    )
assert "accepted_schwarz_state=bounded_schwarz_state" in compact, (
    "accepted coefficients/state must be carried to the next lambda point"
)
assert "bounded_schwarz_state=accepted_schwarz_state" in compact, (
    "failed points must restore the last collectively accepted coefficients/state"
)
continuation_loop = compact.split("dowhile(.not.interface_continuation%finished)", 1)[1]
failure_branch = continuation_loop.split("if(.not.collective_ok)then", 1)[1].split("error stop", 1)[0]
assert "bounded_schwarz_state=accepted_schwarz_state" in failure_branch

diagnostic_call = "call record_dg_hybrid_interface_continuation_diagnostic"
assert route.count(diagnostic_call) == 2, (
    "solver failure and measured point paths must each have one diagnostic/acceptance call"
)
rollback = continuation_loop.split("if(.not.collective_ok)then", 1)[1].split("error stop", 1)[0]
assert diagnostic_call.replace(" ", "") in rollback and ".false.,diagnostic_ok" in rollback, (
    "a rejected point must record rollback before stopping"
)
measured_path = route.split("if(.not.collective_ok)then", 1)[1].split("endif", 1)[1]
assert diagnostic_call in measured_path and ".true." in measured_path, (
    "a solver-successful point must be measured before certification"
)
assert measured_path.index(diagnostic_call) < measured_path.index("accepted_schwarz_state="), (
    "accepted coefficient snapshot must follow diagnostic certification"
)

diagnostic = source.split(
    "subroutine record_dg_hybrid_interface_continuation_diagnostic", 1
)[1].split("end subroutine record_dg_hybrid_interface_continuation_diagnostic", 1)[0]
for field in (
    "lambda=",
    "diagnostic_state_lambda=",
    "accepted_cg_steps=",
    "residual=",
    "orthogonality_defect=",
    "electron_defect=",
    "rayleigh_energy_trace=",
    "scaled_interface_action_norm=",
    "measurement_status=",
    "status=",
    "continuation_fingerprint=",
):
    assert field in diagnostic, f"continuation record is missing {field}"
assert "apply_dg_hybrid_schwarz_hamiltonian(" in diagnostic, (
    "Rayleigh trace must use the live H(lambda) operator"
)
assert "call apply_dg_hybrid_schwarz_h(" not in diagnostic, (
    "diagnostic H action must not mutate the production peer-exchange counter"
)
assert "apply_dg_hybrid_schwarz_s(" in diagnostic, (
    "Rayleigh normalization must use the live metric callback"
)
assert "bounded_interface_scale*bounded_fixed_payload%interface_rows" in re.sub(
    r"\s+", "", diagnostic
), "interface diagnostic must use the scaled complete SIPG row block"
assert "mpi_allreduce" in diagnostic and "mpi_sum" in diagnostic, (
    "row-owned diagnostic contributions must be reduced collectively"
)
assert "continuation_fingerprint_min" in diagnostic and "continuation_fingerprint_max" in diagnostic, (
    "the recorded continuation fingerprint must be rank consistent"
)
assert "measurement_available=collective_diagnostic_ok" in re.sub(r"\s+", "", diagnostic), (
    "operator-measurement failure must select a collective fallback, not suppress the record"
)
assert "rayleigh_energy_trace=huge(1d0)" in re.sub(r"\s+", "", diagnostic)
assert "scaled_interface_action_norm=huge(1d0)" in re.sub(r"\s+", "", diagnostic)
diagnostic_compact = re.sub(r"\s+", "", diagnostic)
assert diagnostic_compact.index("measurement_available=collective_diagnostic_ok") < diagnostic_compact.index(
    "callaccept_dg_hybrid_interface_point"
), "measurement must be known collectively before continuation acceptance"
assert "callvalidate_dg_hybrid_schwarz_dynamic_receipt" in diagnostic_compact
assert diagnostic_compact.index("callvalidate_dg_hybrid_schwarz_dynamic_receipt") < diagnostic_compact.index(
    "callaccept_dg_hybrid_interface_point"
), "common dynamic Schwarz state must be certified before continuation acceptance"
for certification in (
    "callmpi_allreduce(bounded_last_accepted_cg_steps",
    "callmpi_allreduce(continuation_fingerprint",
    "callmpi_allreduce(merge(1,0,local_accept),minimum_accept_request",
    "callmpi_allreduce(merge(1,0,local_accept),maximum_accept_request",
    "callmpi_allreduce(diagnostic_state_lambda,minimum_state_lambda",
    "callmpi_allreduce(diagnostic_state_lambda,maximum_state_lambda",
    "continuation_fingerprint==continuation_state%fingerprint",
    "minimum_steps>=0",
    "maximum_steps<=dg_hybrid_fragment_cg_steps",
    "callcomm_logical_and(record_ok",
):
    assert certification in diagnostic_compact, f"missing pre-accept record certification: {certification}"
    assert diagnostic_compact.index(certification) < diagnostic_compact.index(
        "callaccept_dg_hybrid_interface_point"
    ), f"record certification occurs after continuation mutation: {certification}"
assert "local_accept.and.measurement_available" in diagnostic_compact, (
    "unavailable diagnostics must never reach a true continuation acceptance"
)
for value in (
    "divided_fragment_residual",
    "divided_fragment_orthogonality",
    "bounded_schwarz_state%electron_defect",
):
    assert f"finite_diagnostic_value({value})" in re.sub(r"\s+", "", diagnostic), (
        f"rollback records must retain a finite sentinel for nonfinite {value}"
    )
assert "rank_local==0" in re.sub(r"\s+", "", diagnostic), (
    "only rank zero may print the deterministic continuation record"
)
for forbidden in ("mpi_allgather", "mpi_allgatherv"):
    assert forbidden not in diagnostic, "diagnostics must not gather coefficient matrices"

h_callback = source.split("subroutine apply_dg_hybrid_schwarz_h(", 1)[1].split(
    "end subroutine apply_dg_hybrid_schwarz_h", 1
)[0]
assert "bounded_interface_scale" in h_callback, (
    "Schwarz Hamiltonian callback must pass the live lambda"
)
preconditioner = source.split(
    "subroutine assemble_dg_hybrid_schwarz_local_preconditioner_blocks", 1
)[1].split("end subroutine assemble_dg_hybrid_schwarz_local_preconditioner_blocks", 1)[0]
assert re.search(
    r"bounded_interface_scale\s*\*\s*bounded_fixed_payload%interface_rows",
    preconditioner,
), "local preconditioner must scale the complete interface rows by the same lambda"

print("fixed-DC-density DG interface continuation route: PASS")
