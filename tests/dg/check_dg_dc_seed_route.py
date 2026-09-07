#!/usr/bin/env python3
"""Source contract for the conventional-DC seed route."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MAIN = (ROOT / "src/gs/main_dft.f90").read_text()
CHECKPOINT = (ROOT / "src/gs/dc/dg_dc_seed_checkpoint.f90").read_text()
DCDFT = (ROOT / "src/gs/dc/dcdft.f90").read_text()
SCF = (ROOT / "src/gs/scf_iteration_dft.f90").read_text()
INITIALIZATION_DFT = (ROOT / "src/gs/initialization_dft.f90").read_text()
MAIN_DFT_MD = (ROOT / "src/gs/main_dft_md.f90").read_text()
CANONICAL_PP = (ROOT / "src/gs/dc/dg_canonical_pp_fingerprint.f90").read_text()


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text.lower())


def call_extent(source: str, name: str, start: int = 0) -> tuple[int, int, str]:
    match = re.search(rf"\bcall\s+{name}\s*\(", source[start:], re.I)
    assert match, f"missing production call to {name}"
    begin = start + match.start()
    cursor = start + match.end()
    depth = 1
    while cursor < len(source) and depth:
        if source[cursor] == "(":
            depth += 1
        elif source[cursor] == ")":
            depth -= 1
        cursor += 1
    assert depth == 0, f"unterminated production call to {name}"
    return begin, cursor, source[begin:cursor]


def subroutine_extent(source: str, name: str) -> tuple[int, int, str]:
    match = re.search(rf"\bsubroutine\s+{name}\s*\(", source, re.I)
    assert match, f"missing production subroutine {name}"
    ending = re.search(rf"\bend\s+subroutine\s+{name}\b", source[match.end():], re.I)
    assert ending, f"unterminated production subroutine {name}"
    end = match.end() + ending.end()
    return match.start(), end, source[match.start():end]


def require_collective_ok_gate(source: str, begin: int, end: int, label: str) -> None:
    region = compact(source[begin:end])
    assert (
        "callcomm_logical_and(dg_dc_seed_ok,dg_dc_seed_collective_ok,dc%icomm_tot)"
        in region
    ), f"{label} local success is not reduced across the total communicator"
    assert re.search(
        r"if\(\.not\.dg_dc_seed_collective_ok\).*?errorstop", region, re.I
    ), f"{label} failure does not terminate collectively before the next collective"


checkpoint_specification = CHECKPOINT[:CHECKPOINT.lower().index("contains")].replace("&\n", "")
for public_name in (
    "resolve_dg_dc_seed_mode",
    "build_dg_dc_seed_contract",
    "restore_dg_dc_seed_payload",
):
    assert re.search(rf"\bpublic\s*::[^\n]*\b{public_name}\b", checkpoint_specification, re.I), (
        f"DC seed state helper is not public: {public_name}"
    )

early_band_return_position = MAIN.lower().index("if(theory=='dft_band'.and.iperiodic/=3) return")
preinit_gate_position = MAIN.lower().rfind(
    "if(trim(dg_dc_seed_mode)/='off')", 0, early_band_return_position
)
assert preinit_gate_position >= 0, (
    "enabled seed scope must be rejected before the nonperiodic-band early return and DC init"
)
preinit_gate = compact(MAIN[preinit_gate_position:early_band_return_position])
for globally_required_scope in (
    "yn_dc",
    "yn_dg_dc_overlapping_wannier",
    "theory",
    "yn_spinorbit",
    "plus_u_on",
    "yn_hse",
    "yn_fix_func",
    "yn_jm",
):
    assert globally_required_scope in preinit_gate, (
        f"pre-init seed scope gate omits {globally_required_scope}"
    )
assert (
    "callcomm_logical_and(dg_dc_seed_ok,dg_dc_seed_collective_ok,nproc_group_global)"
    in preinit_gate
), "pre-init seed scope is not agreed on the initialized global communicator"
assert re.search(r"if\(\.not\.dg_dc_seed_collective_ok\).*?errorstop", preinit_gate), (
    "unsupported pre-init seed scope is not rejected on every rank"
)

initialization_position = MAIN.lower().index("call initialization2_dft")
scf_position = MAIN.lower().index("call scf_iteration_dft")
probe_position, _, probe_call = call_extent(MAIN, "probe_dg_dc_seed")
assert initialization_position < probe_position < scf_position, (
    "seed probing must occur after initialization2_dft and before conventional SCF"
)
for exact_total_topology in ("dc%icomm_tot", "dc%id_tot", "dc%isize_tot"):
    assert exact_total_topology in MAIN[initialization_position:scf_position].lower(), (
        f"seed route must use the preserved total-system topology: {exact_total_topology}"
    )
assert "dg_dc_seed_directory" in probe_call.lower(), (
    "seed probe does not use the user-selected seed directory"
)

seed_gate_begin = MAIN.lower().index("if(trim(dg_dc_seed_mode)/='off')", initialization_position)
seed_gate_end = MAIN.lower().index("call prepare_dg_dc_seed_contract_inputs", seed_gate_begin)
seed_gate = compact(MAIN[seed_gate_begin:seed_gate_end])
for unsupported_hamiltonian in ("plus_u_on", "yn_hse", "yn_fix_func", "yn_jm"):
    assert unsupported_hamiltonian in seed_gate, (
        "the strict seed entry gate does not reject unsupported Hamiltonian control: "
        f"{unsupported_hamiltonian}"
    )
require_collective_ok_gate(
    MAIN, seed_gate_begin, seed_gate_end, "supported-scope entry validation"
)

mode_region = MAIN[initialization_position:scf_position]
assert re.search(r"select\s+case\s*\(\s*trim\s*\(\s*dg_dc_seed_mode\s*\)\s*\)", mode_region, re.I), (
    "missing explicit off/write/read/auto DC seed state machine"
)
for mode in ("off", "write", "read", "auto"):
    assert re.search(rf"case\s*\(\s*['\"]{mode}['\"]\s*\)", mode_region, re.I), (
        f"DC seed state machine omits mode={mode}"
    )
assert re.search(r"call\s+resolve_dg_dc_seed_mode\b", mode_region, re.I), (
    "production state machine must share the tested mode resolver"
)

contract_gate_position, contract_gate_end, _ = call_extent(
    MAIN, "build_dg_dc_seed_contract"
)
state_machine_position = MAIN.lower().index(
    "select case(trim(dg_dc_seed_mode))", contract_gate_end
)
require_collective_ok_gate(
    MAIN, contract_gate_end, state_machine_position, "preserved total topology validation"
)
topology_gate = compact(MAIN[contract_gate_end:state_machine_position])
for exact_topology_match in (
    "dc%id_tot==dg_dc_seed_contract%rank",
    "dc%isize_tot==dg_dc_seed_contract%mpi_size",
):
    assert exact_topology_match in topology_gate, (
        f"preserved total topology does not prove {exact_topology_match}"
    )

scf_guard = MAIN[max(0, scf_position - 300):scf_position]
assert re.search(r"if\s*\(\s*dg_dc_seed_run_scf\s*\)\s*then", scf_guard, re.I), (
    "valid read/auto seeds do not guard and skip conventional scf_iteration_dft"
)
assert "DC #SCF =" in (ROOT / "src/gs/scf_iteration_dft.f90").read_text(), (
    "test contract lost the conventional DC SCF receipt"
)
_, _, scf_call = call_extent(MAIN, "scf_iteration_dft")
for seed_call_argument in (
    "dg_dc_seed_publish",
    "dg_dc_seed_electron_tolerance",
):
    assert seed_call_argument in scf_call.lower(), (
        f"the conventional SCF seed convergence call omits {seed_call_argument}"
    )
for non_seed_source, label in (
    (INITIALIZATION_DFT, "initial MD SCF"),
    (MAIN_DFT_MD, "MD-step SCF"),
):
    _, _, non_seed_call = call_extent(non_seed_source, "scf_iteration_dft")
    non_seed_compact = compact(non_seed_call)
    assert ",.false.,0d0" in non_seed_compact, (
        f"{label} does not explicitly disable the strict DG DC seed convergence gate"
    )
_, _, scf_body = subroutine_extent(SCF, "scf_iteration_dft")
scf_compact = compact(scf_body)
for seed_charge_gate in (
    "require_dg_dc_seed_electron_count",
    "required_dg_dc_seed_electron_count_tolerance",
    "sum(dc%rho_tot_s(1)%f)",
    "dc%system_tot%hvol",
    "dc%elec_num_tot",
    "call comm_summation",
):
    assert compact(seed_charge_gate) in scf_compact, (
        f"strict seed SCF convergence omits {seed_charge_gate}"
    )
assert "[DG-DC-SEED-WAIT]" in scf_body, (
    "continued SCF for seed electron-count convergence lacks an audit receipt"
)
assert re.search(
    r"if\s*\(\s*sum1\s*<\s*threshold\s*\.and\.\s*"
    r"dg_dc_seed_electron_count_converged\s*\)\s*then",
    scf_body,
    re.I,
), "ordinary residual convergence can still bypass strict seed electron-count convergence"

read_call_position, _, _ = call_extent(MAIN, "read_dg_dc_seed")
restore_call_position, restore_call_end, restore_call = call_extent(
    MAIN, "restore_dg_dc_seed_payload"
)
validate_call_position, _, _ = call_extent(MAIN, "validate_dg_dc_seed_state")
rebuild_call_position, _, _ = call_extent(MAIN, "rebuild_dg_dc_seed_derived_state_dcdft")
assert (
    probe_position
    < read_call_position
    < restore_call_position
    < validate_call_position
    < rebuild_call_position
    < scf_position
), (
    "valid seed state must be read and restored before the guarded conventional SCF"
)
for exact_state in (
    "spsi%rwf",
    "dc%rho_tot_s(1)%f",
    "dc%vloc_tot(1)%f",
    "energy%esp",
    "system%rocc",
    "system%mu",
    "sum1",
    "miter",
):
    assert exact_state in restore_call.lower(), f"production restore omits {exact_state}"
require_collective_ok_gate(
    MAIN, restore_call_end, validate_call_position, "restored-payload allocation"
)
assert re.search(r"if\s*\(\s*dg_dc_seed_fatal\s*\).*?error\s+stop", mode_region, re.I | re.S), (
    "read-absent and read/auto-invalid seed states must terminate collectively"
)
state_machine_fatal_position = mode_region.lower().index("if(dg_dc_seed_fatal)")
state_machine_load_position = mode_region.lower().index(
    "if(dg_dc_seed_load)", state_machine_fatal_position
)
state_machine_failure_region = mode_region[
    state_machine_fatal_position:state_machine_load_position
]
assert "dg_dc_seed_message" in state_machine_failure_region.lower(), (
    "seed state-machine failure discards the probe/resolver diagnostic"
)
assert "[DG-DC-SEED-ERROR]" in state_machine_failure_region, (
    "seed state-machine failure lacks a recognizable diagnostic receipt"
)
assert re.search(
    r"if\s*\(\s*dc%id_tot\s*==\s*0\s*\).*?write\s*\(\s*error_unit",
    state_machine_failure_region,
    re.I | re.S,
), "seed state-machine diagnostic is not emitted by total-system rank zero"
assert re.search(r"flush\s*\(\s*error_unit\s*\)", state_machine_failure_region, re.I), (
    "seed state-machine diagnostic is not flushed before fatal termination"
)
assert (
    state_machine_failure_region.index("[DG-DC-SEED-ERROR]")
    < state_machine_failure_region.lower().index("error stop")
), "seed state-machine terminates before publishing its detailed diagnostic"

convergence_position = MAIN.lower().index("if(.not.(sum1<threshold))")
write_position, write_call_end, write_call = call_extent(MAIN, "write_dg_dc_seed")
dispatch_positions = [
    position
    for token in (
        "call run_dg_hybrid_continuation_ground_state_for_main",
        "call run_dg_overlapping_wannier_ground_state_for_main",
        "call dc_lcfo_flux",
        "call dc_lcfo(",
    )
    if (position := MAIN.lower().find(token, convergence_position)) >= 0
]
assert dispatch_positions, "cannot locate any post-DC local-basis dispatch"
assert convergence_position < write_position < min(dispatch_positions), (
    "DC seed must be published only after convergence and before OW/Hybrid/LCFO dispatch"
)
publish_guard_position = MAIN.lower().rfind(
    "if(dg_dc_seed_publish)then", convergence_position, write_position
)
assert publish_guard_position >= convergence_position, (
    "seed write is not guarded by the write/auto publication decision"
)
capture_position, capture_call_end, capture_call = call_extent(
    MAIN, "capture_dg_dc_seed_payload_dcdft"
)
assert convergence_position < capture_position < write_position, (
    "the converged production state must be captured immediately before publication"
)
for exact_state in ("system", "energy", "spsi", "dc", "sum1", "miter"):
    assert exact_state in capture_call.lower(), f"production capture omits {exact_state}"
require_collective_ok_gate(
    MAIN, capture_call_end, write_position, "captured-payload allocation"
)
publish_failure_region = MAIN[write_call_end:write_call_end + 700]
assert "dg_dc_seed_message" in publish_failure_region.lower(), (
    "production seed publication failure discards the checkpoint writer's diagnostic"
)
assert "[DG-DC-SEED-ERROR]" in publish_failure_region, (
    "production seed publication failure lacks a recognizable diagnostic receipt"
)
assert re.search(
    r"if\s*\(\s*\.not\.\s*dg_dc_seed_ok\s*\).*?error\s+stop",
    publish_failure_region,
    re.I | re.S,
), "production seed publication failure is not fatal after reporting its diagnostic"

_, _, capture_body = subroutine_extent(DCDFT, "capture_dg_dc_seed_payload_dcdft")
for exact_copy in (
    "payload%rwf=spsi%rwf",
    "payload%rho_tot=dc%rho_tot_s(1)%f",
    "payload%vloc_tot=dc%vloc_tot(1)%f",
    "payload%esp=energy%esp",
    "payload%rocc=system%rocc",
    "payload%mu=system%mu",
    "payload%residual=residual",
    "payload%iteration=iteration",
):
    assert compact(exact_copy) in compact(capture_body), f"production capture omits {exact_copy}"

_, _, rebuild_body = subroutine_extent(DCDFT, "rebuild_dg_dc_seed_derived_state_dcdft")
for rebuild_step in (
    "dc%rho_tot%f=dc%rho_tot%f+dc%rho_tot_s(ispin)%f",
    "call calc_density",
    "rho%f=rho%f+rho_s(ispin)%f",
    "call calc_vlocal_fragment_dcdft",
):
    assert compact(rebuild_step) in compact(rebuild_body), (
        f"restored-state rebuild omits {rebuild_step}"
    )
for forbidden_overwrite in (
    "call ne2mu_dcdft",
    "call calc_rho_total_dcdft",
    "call update_density_and_potential",
    "call copy_density",
):
    assert compact(forbidden_overwrite) not in compact(rebuild_body), (
        f"restored-state rebuild overwrites the exact seed via {forbidden_overwrite}"
    )

assert MAIN.count("[DG-DC-SEED]") == 1, (
    "the route must emit exactly one canonical DC seed receipt"
)
receipt_position = MAIN.index("[DG-DC-SEED]")
receipt = MAIN[receipt_position:receipt_position + 500].lower()
for field in ("mode=", "publication_id=", "scf_skipped=", "mpi_size=", "mapping_fingerprint="):
    assert field in receipt, f"DC seed receipt omits {field}"

prepare_position, prepare_end, _ = call_extent(MAIN, "prepare_dg_dc_seed_contract_inputs")
contract_position, _, contract_call = call_extent(MAIN, "build_dg_dc_seed_contract")
require_collective_ok_gate(
    MAIN, prepare_end, contract_position, "production contract preparation"
)
contract_call = compact(contract_call)
for downstream_control in (
    "dg_ow_",
    "wannier_",
    "wannierpw",
    "lcfo",
    "window",
    "yn_dg_hybrid",
    "yn_dg_dc_overlapping_wannier",
    "dg_hybrid_divided_mixing",
):
    assert downstream_control not in contract_call, (
        "downstream localization/W90/PW/LCFO/window controls must not enter the "
        f"conventional DC seed fingerprint: {downstream_control}"
    )

_, _, convergence_fingerprint_body = subroutine_extent(
    MAIN.replace(
        "integer(int64) function dg_dc_seed_convergence_fingerprint()result(hash)",
        "subroutine dg_dc_seed_convergence_fingerprint()",
        1,
    ).replace(
        "end function dg_dc_seed_convergence_fingerprint",
        "end subroutine dg_dc_seed_convergence_fingerprint",
        1,
    ),
    "dg_dc_seed_convergence_fingerprint",
)
convergence_fingerprint_compact = compact(convergence_fingerprint_body)
assert "hash_character(hash,trim(method_mixing))" in convergence_fingerprint_compact
assert "hash_real(hash,mixing%mixrate)" in convergence_fingerprint_compact
assert "dg_hybrid_divided_mixing" not in convergence_fingerprint_compact, (
    "post-seed Hybrid mixing selector changed conventional DC seed compatibility"
)

_, _, operator_fingerprint_body = subroutine_extent(
    MAIN.replace(
        "integer(int64) function dg_dc_seed_operator_input_fingerprint()result(hash)",
        "subroutine dg_dc_seed_operator_input_fingerprint()",
        1,
    ).replace(
        "end function dg_dc_seed_operator_input_fingerprint",
        "end subroutine dg_dc_seed_operator_input_fingerprint",
        1,
    ),
    "dg_dc_seed_operator_input_fingerprint",
)
operator_fingerprint_compact = compact(operator_fingerprint_body)
assert "canonical_pp_fingerprint(pp)" in operator_fingerprint_compact, (
    "the production DC seed does not use the canonical pseudopotential fingerprint"
)
assert "if(pp_fingerprint==0_int64)thenhash=0_int64;return" in operator_fingerprint_compact, (
    "an incomplete canonical pseudopotential can enter the seed contract"
)
assert "hash_pp_info_for_dg_dc_seed" not in MAIN.lower(), (
    "the seed route still exposes the old whole-allocation pseudopotential hash"
)
for canonical_contract in (
    "SALMON-DG-CANONICAL-PP",
    "do radial=1,pp%nrloc(element)",
    "do radial=1,pp%nrps(element)",
    "do radial=1,pp%nrps_ao(element)",
    "nlcc_meaning_end(pp,element)",
):
    assert compact(canonical_contract) in compact(CANONICAL_PP), (
        f"canonical pseudopotential contract omits {canonical_contract}"
    )
for forbidden_scratch in (
    "pp%zion",
    "pp%vpp=",
    "pp%upp=",
    "pp%dvpp=",
    "pp%dupp=",
    "pp%vpp_f",
    "pp%upp_f",
    "pp%dupptbl_ao",
):
    assert forbidden_scratch not in CANONICAL_PP.lower(), (
        f"canonical pseudopotential hash consumes reader scratch: {forbidden_scratch}"
    )
assert "canonical_pp_valence_sum(pp)" in compact(MAIN), (
    "the Hybrid checkpoint pseudopotential receipt still consumes scalar pp%zion scratch"
)

assert re.search(
    r"ownership_map\s*\(\s*8\s*\)\s*=\s*int\s*\(\s*dc%id_tot\s*\+\s*1\s*,\s*int64\s*\)",
    MAIN,
    re.I,
), "zero-based total rank must be encoded nonzero before the contract zero-sentinel check"

print("PASS conventional DC seed production route contract")
