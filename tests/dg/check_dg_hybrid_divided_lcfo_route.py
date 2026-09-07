#!/usr/bin/env python3
"""Production contract for fixed-density continuation and one terminal LCFO solve."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()

ENTRY_NAME = "subroutine run_dg_hybrid_divided_ground_state_for_main"
ENTRY_END = "end subroutine run_dg_hybrid_divided_ground_state_for_main"
assert ENTRY_NAME in SOURCE and ENTRY_END in SOURCE
entry = SOURCE[SOURCE.index(ENTRY_NAME) : SOURCE.index(ENTRY_END)]

assert entry.count("call initialize_dg_hybrid_interface_continuation") == 1, (
    "production route must initialize exactly one fixed-density interface continuation"
)
assert "do while(.not.interface_continuation%finished)" in entry
assert entry.count("call accept_dg_hybrid_interface_point") == 0, (
    "the driver must not advance continuation before diagnostic certification"
)
scf_position = entry.index("call initialize_dg_hybrid_interface_continuation")
solve_position = entry.find("call solve_dg_hybrid_generalized_once_and_publish")
assert solve_position > scf_position, "one terminal LCFO solve must follow interface continuation"
assert entry.count("call solve_dg_hybrid_generalized_once_and_publish") == 1
guard_position = entry.find("call validate_dg_hybrid_schwarz_dynamic_receipt")
row_position = entry.find("final_hrows=bounded_fixed_payload%kinetic_rows")
assert scf_position < guard_position < row_position < solve_position, (
    "terminal continuation guard must precede full row composition and the one LCFO solve"
)
terminal_guard = entry[guard_position:row_position]
for required in (
    "interface_continuation%finished",
    "interface_continuation%lambda==1d0",
    "bounded_interface_scale==1d0",
    "accepted_interface_scale==1d0",
    "interface_continuation%fingerprint/=0_int64",
    "bounded_schwarz_state%fingerprint/=0_int64",
    "interface_continuation%basis_generation==bounded_schwarz_state%basis_generation",
    "interface_continuation%mapping_fingerprint==bounded_schwarz_state%mapping_fingerprint",
    "call validate_dg_hybrid_schwarz_dynamic_receipt",
):
    assert required in terminal_guard, f"terminal guard is missing: {required}"
assert "call comm_logical_and" in terminal_guard
assert "mpi_allreduce" in terminal_guard, (
    "terminal fingerprints must be checked for rank consistency"
)
solve_call = entry[solve_position:].split("ok,message)", 1)[0]
assert "call mpi_allreduce(frame_fingerprint,global_frame_fingerprint" in entry, (
    "fragment-local reference-frame receipts must be reduced to one global receipt"
)
assert "final_operator_fingerprint,global_frame_fingerprint" in solve_call, (
    "terminal state validation must receive the collective reference-frame receipt"
)
assert "electronic_temperature=max(0d0,temperature)" in solve_call, (
    "terminal occupations must receive SALMON's atomic-unit electronic temperature"
)
assert "electronic_temperature=bounded_schwarz_state%temperature" not in solve_call, (
    "Schwarz state temperature is in kelvin and must not enter the hartree occupation kernel"
)
assert "occupation_electron_tolerance=dg_dc_gs_electron_count_tolerance" in solve_call

terminal = entry[scf_position:]
potential_position = terminal.find("call extract_dg_hybrid_core_local_potential")
projection_position = terminal.find("call assemble_dg_hybrid_local_potential_rows")
fingerprint_position = terminal.find("call ow_fingerprint_distributed_matrix")
assert 0 < potential_position < projection_position < fingerprint_position, (
    "terminal LCFO must project the once-refreshed converged potential"
)
assert terminal.count("call assemble_dg_hybrid_local_potential_rows") == 1, (
    "terminal density-independent and converged-potential rows must be composed once"
)
assert "bounded_fixed_payload%kinetic_rows+bounded_fixed_payload%nonlocal_rows+&" in terminal
assert "bounded_fixed_payload%interface_rows+final_local_potential_rows" in terminal
assert "final_srows=bounded_fixed_payload%metric_rows" in terminal
assert "fixed-density/non-self-consistent" in entry[solve_position:], (
    "terminal output must label the LCFO result as fixed-density/non-self-consistent"
)

checkpoint_position = entry.find("call write_rt_dg_hybrid_occupied_checkpoint")
assert checkpoint_position > solve_position, "occupied checkpoint must follow terminal LCFO"
assert entry.count("call write_rt_dg_hybrid_occupied_checkpoint") == 1
post_lcfo = entry[solve_position:]
for forbidden in (
    "call run_dg_hybrid_divided_scf",
    "call mix_dg_hybrid_divided_density",
    "call update_dg_hybrid_divided_potential",
    "call dg_dc_update_potential_from_distributed_density",
    "reconstruct_dg_hybrid_density",
    "post_lcfo_density",
):
    assert forbidden not in post_lcfo, f"post-LCFO density work is forbidden: {forbidden}"

SOLVER_NAME = "subroutine solve_dg_hybrid_schwarz_fragments"
SOLVER_END = "end subroutine solve_dg_hybrid_schwarz_fragments"
assert SOLVER_NAME in SOURCE and SOLVER_END in SOURCE
solver = SOURCE[SOURCE.index(SOLVER_NAME) : SOURCE.index(SOLVER_END)]
assert "solve_dg_hybrid_generalized" not in solver, (
    "fragment continuation points must not perform a full generalized eigensolve"
)

for forbidden in (
    "allocate(ow_hybrid_hrows(ntarget,ntarget)",
    "allocate(ow_hybrid_srows(ntarget,ntarget)",
    "complex(8)::hybrid_h(ntarget,ntarget)",
    "complex(8)::hybrid_s(ntarget,ntarget)",
):
    assert forbidden not in entry.replace(" ", ""), (
        f"divided route must retain row-distributed LCFO storage: {forbidden}"
    )

print("fixed-density divided WF+PW LCFO route contract: PASS")
