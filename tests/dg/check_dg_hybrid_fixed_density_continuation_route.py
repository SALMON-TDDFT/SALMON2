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
assert re.search(r"call\s+accept_dg_hybrid_interface_point\s*\(", route), (
    "successful points must advance the collective continuation state"
)

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
