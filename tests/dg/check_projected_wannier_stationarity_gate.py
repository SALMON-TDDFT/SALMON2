#!/usr/bin/env python3
"""Protect divided-GS v5 publication and current Hybrid RT stationarity gates."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RT_SOURCE = (ROOT / "src/rt/main_tddft.f90").read_text(errors="replace").lower()
GS_SOURCE = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()


def subroutine(source: str, name: str) -> str:
    start = source.find(f"subroutine {name}")
    if start < 0:
        return ""
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
    hybrid_gs = subroutine(gs_source, "publish_dg_hybrid_divided_v5")
    if not ordered(
        hybrid_gs,
        (
            "call collective_rt_dg_hybrid_publication_precondition",
            "terminal divided v5 publication precondition failed:",
            "call collective_rt_dg_hybrid_publication_mapping_precondition",
            "terminal divided v5 row mapping failed:",
            "call validate_dg_hybrid_v5_publication_rank_policy",
            "authorization%valid=publication_authorized",
            "call publish_rt_dg_hybrid_checkpoint_v5",
            "terminal divided v5 write failed:",
        ),
    ):
        problems.append("divided GS must validate and authorize before checkpoint publication")
    return problems


assert not violations(RT_SOURCE, GS_SOURCE), violations(RT_SOURCE, GS_SOURCE)

# Mutation fixtures cover acceptance and publication ordering, not just token
# presence.  Each mutation must invalidate the corresponding safety contract.
for token in (
    "call initialize_rt_dg_hybrid_from_checkpoint",
    "error stop 'hybrid dg rt initialization failed'",
    "call initialize_rt_dg_hybrid_stationarity",
    "call propagate_rt_dg_hybrid_length_gauge",
    "call evaluate_rt_dg_hybrid_stationarity",
    "error stop 'hybrid dg rt zero-field stationarity failed'",
):
    mutation = RT_SOURCE.replace(token, "omitted", 1)
    assert violations(mutation, GS_SOURCE), token
for token in (
    "call collective_rt_dg_hybrid_publication_precondition",
    "terminal divided v5 publication precondition failed:",
    "call collective_rt_dg_hybrid_publication_mapping_precondition",
    "terminal divided v5 row mapping failed:",
    "call validate_dg_hybrid_v5_publication_rank_policy",
    "authorization%valid=publication_authorized",
    "call publish_rt_dg_hybrid_checkpoint_v5",
    "terminal divided v5 write failed:",
):
    mutation = GS_SOURCE.replace(token, "omitted", 1)
    assert violations(RT_SOURCE, mutation), token

print("PASS current v5 Hybrid RT stationarity and GS publication gates with mutation fixtures")
