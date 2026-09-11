#!/usr/bin/env python3
"""Static reachability contract for the formal Hybrid continuation -> v5 route."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/main_dft.f90").read_text()
CHECKPOINT = (ROOT / "src/rt/dg/rt_dg_hybrid_checkpoint_v5.f90").read_text()


def compact(text: str) -> str:
    return "".join(text.lower().split()).replace("&", "")


def routine(text: str, name: str) -> str:
    lower = text.lower()
    start = lower.index(f"subroutine {name}")
    end = lower.index("end subroutine", start)
    return text[start:end]


route = routine(SOURCE, "run_dg_hybrid_continuation_ground_state_for_main")
continuation = routine(SOURCE, "run_dg_hybrid_divided_ground_state_for_main")
publisher = routine(SOURCE, "publish_dg_hybrid_divided_v5")
route_c = compact(route)
continuation_c = compact(continuation)
publisher_c = compact(publisher)

def require_publication_rank_policy(body: str) -> None:
    assert body.count("callvalidate_dg_hybrid_v5_publication_rank_policy") == 1
    assert "if(present(publication_receipt))rank_policy_receipt=publication_receipt" in body
    assert "dg_hybrid_symmetry_energy_window,requested_rank,certified_rank,n" in body

require_publication_rank_policy(publisher_c)
for old in (
    "callvalidate_dg_hybrid_v5_publication_rank_policy",
    "if(present(publication_receipt))rank_policy_receipt=publication_receipt",
):
    mutated = publisher_c.replace(old, "removed_rank_policy", 1)
    try:
        require_publication_rank_policy(mutated)
    except AssertionError:
        pass
    else:
        raise AssertionError(f"publication-rank policy mutation survived: {old}")

# The compatibility selector enters the same local-plus-terminal production driver.
assert "callrun_dg_hybrid_divided_ground_state_for_main" in route_c
assert "run_dg_hybrid_concrete_continuation" not in SOURCE.lower()
assert continuation_c.count("callpublish_dg_hybrid_divided_v5(") == 1
publish = continuation_c.index("callpublish_dg_hybrid_divided_v5(")
for token in ("enddoterminal_lcfo_refinement", "callwrite_rt_dg_hybrid_occupied_checkpoint",
              "callobserve_dg_hybrid_terminal_refinement"):
    assert continuation_c.index(token) < publish
assert "terminal_refinement=terminal_refinement_receipt" in continuation_c[publish:]
assert "callsolve_dg_hybrid_generalized_once_and_publish" not in continuation_c[publish:]
assert "calldg_dc_update_potential_from_distributed_density" not in continuation_c[publish:]

# v5 preserves localized construction rows and occupied LCFO amplitudes; it
# never rotates them into replicated spectral basis state.
for token in (
    "publish_rt_dg_hybrid_checkpoint_v5",
    "payload%initial_occupied_amplitudes",
    "payload%basis_point_offsets",
    "payload%basis_support_ids",
    "payload%operator_offsets",
    "payload%operator_columns",
    "dc%ppg_tot%nlma",
    "payload%energy_receipt=[checkpoint_energy%e_tot",
    "callcalc_total_energy_periodic",
):
    assert token in publisher_c, token
for forbidden in (
    "full_coefficients",
    "full_metric",
    "collect_dg_hybrid_full_rows",
    "write_rt_dg_hybrid_ground_state_checkpoint",
    "publish_dg_hybrid_divided_v3",
):
    assert forbidden not in continuation_c
    assert forbidden not in publisher_c

# The remaining checkpoint component is publication preflight plus occupied
# v2 compatibility only; no dense v3 payload/writer/reader implementation.
checkpoint_c = compact(CHECKPOINT)
for forbidden in (
    "s_rt_dg_hybrid_ground_state_payload",
    "write_rt_dg_hybrid_ground_state_checkpoint",
    "read_rt_dg_hybrid_ground_state_checkpoint",
    "full_coefficients",
    "full_metric",
):
    assert forbidden not in checkpoint_c

# Removing the reachable publication must invalidate the publication count.
mutated = continuation_c.replace("callpublish_dg_hybrid_divided_v5(", "callremoved_v5(", 1)
assert mutated.count("callpublish_dg_hybrid_divided_v5(") == 0
print("PASS Hybrid compatibility selector reaches terminal distributed v5 publication")
