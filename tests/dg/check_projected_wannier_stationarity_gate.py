#!/usr/bin/env python3
"""Protect current checkpoint and coefficient-space RT stationarity gates.

The prior inventory test named an obsolete projected-Wannier seed receipt.  The
bridge already used certified hybrid checkpoints and coefficient-space RT.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RT_SOURCE = (ROOT / "src/rt/main_tddft.f90").read_text(errors="replace").lower()
GS_SOURCE = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()


def subroutine(source: str, name: str) -> str:
    start = source.index(f"subroutine {name}")
    end = source.index("end subroutine", start)
    return source[start:end]


def ordered(body: str, tokens: tuple[str, ...]) -> bool:
    positions = [body.find(token) for token in tokens]
    return all(position >= 0 for position in positions) and positions == sorted(positions)


def violations(rt_source: str, gs_source: str) -> list[str]:
    problems: list[str] = []
    hybrid_rt = subroutine(rt_source, "run_dg_hybrid_continuation_rt")
    if not ordered(
        hybrid_rt,
        (
            "call initialize_rt_dg_hybrid_from_checkpoint",
            "error stop 'hybrid dg rt initialization failed'",
            "call initialize_rt_dg_hybrid_stationarity",
            "call propagate_rt_dg_hybrid_length_gauge",
            "call evaluate_rt_dg_hybrid_stationarity",
            "error stop 'hybrid dg rt zero-field stationarity failed'",
        ),
    ):
        problems.append("hybrid checkpoint rejection, propagation, and stationarity gates are out of order")
    coefficient_rt = subroutine(rt_source, "run_dg_overlapping_wannier_coefficient_rt")
    compact_coefficient = "".join(coefficient_rt.replace("&", "").split())
    if "if(.not.ok.or..not.reusable)then" not in compact_coefficient:
        problems.append("coefficient RT must reject a checkpoint that is invalid or not reusable")
    if not ordered(
        coefficient_rt,
        (
            "call read_dg_overlapping_wannier_checkpoint",
            "error stop 'accepted v3 overlapping-wannier checkpoint is required'",
            "call initialize_dg_overlapping_wannier_rt",
            "call advance_dg_overlapping_wannier_rt",
            "call write_dg_overlapping_wannier_rt_restart",
        ),
    ):
        problems.append("coefficient checkpoint acceptance, propagation, and publication are out of order")
    hybrid_gs_start = gs_source.index("if(yn_dg_hybrid_scf=='y')then")
    hybrid_gs_end = gs_source.index("return", hybrid_gs_start)
    hybrid_gs = gs_source[hybrid_gs_start:hybrid_gs_end]
    if not ordered(
        hybrid_gs,
        (
            "call run_dg_hybrid_self_consistent_ground_state",
            "call validate_dg_hybrid_ground_state",
            "ow_hybrid_ground_state%converged=.true.",
            "call write_rt_dg_hybrid_occupied_checkpoint",
            "error stop 'distributed hybrid checkpoint failed'",
        ),
    ):
        problems.append("hybrid GS must validate and mark convergence before checkpoint publication")
    return problems


assert not violations(RT_SOURCE, GS_SOURCE), violations(RT_SOURCE, GS_SOURCE)

# Mutation fixtures cover acceptance and publication ordering, not just token
# presence.  Each mutation must invalidate the corresponding safety contract.
accepts_stale = RT_SOURCE.replace(".not.ok.or..not.reusable", ".not.ok", 1)
assert any("invalid or not reusable" in item for item in violations(accepts_stale, GS_SOURCE))
no_stationarity = RT_SOURCE.replace("call evaluate_rt_dg_hybrid_stationarity", "call omitted_stationarity", 1)
assert any("stationarity gates are out of order" in item for item in violations(no_stationarity, GS_SOURCE))
publishes_unconverged = GS_SOURCE.replace("ow_hybrid_ground_state%converged=.true.", "! mutation removed convergence", 1)
assert any("before checkpoint publication" in item for item in violations(RT_SOURCE, publishes_unconverged))

print("PASS current hybrid/coefficient RT stationarity and checkpoint gates with mutation fixtures")
