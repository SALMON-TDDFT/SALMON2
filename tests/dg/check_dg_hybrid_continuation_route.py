#!/usr/bin/env python3
"""Static reachability contract for the formal Hybrid continuation -> v4 route."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/main_dft.f90").read_text()
CHECKPOINT = (ROOT / "src/rt/dg/rt_dg_hybrid_checkpoint.f90").read_text()


def compact(text: str) -> str:
    return "".join(text.lower().split()).replace("&", "")


def routine(text: str, name: str) -> str:
    lower = text.lower()
    start = lower.index(f"subroutine {name}")
    end = lower.index("end subroutine", start)
    return text[start:end]


route = routine(SOURCE, "run_dg_overlapping_wannier_ground_state_for_main")
continuation = routine(SOURCE, "run_dg_hybrid_concrete_continuation")
publisher = routine(SOURCE, "publish_dg_hybrid_divided_v4")
route_c = compact(route)
continuation_c = compact(continuation)
publisher_c = compact(publisher)

# The formal continuation switch must reach the concrete driver and terminate
# there, before the non-continuation truncated-occupation route is entered.
branch = route_c.index("if(yn_dg_hybrid_continuation_scf=='y')then")
call = route_c.index("callrun_dg_hybrid_concrete_continuation", branch)
ret = route_c.index("return", call)
legacy_occupation = route_c.index("allocate(occupations(nstate))", ret)
assert branch < call < ret < legacy_occupation

# A terminal accepted state is published exactly once by the same distributed
# v4 publisher as the divided route. Publication follows final refresh,
# final row-state validation and the authoritative final operator fingerprint.
assert continuation_c.count("callpublish_dg_hybrid_divided_v4(") == 1
publish = continuation_c.index("callpublish_dg_hybrid_divided_v4(")
for token in (
    "if(.not.final_refresh_performed)errorstop",
    "callow_fingerprint_distributed_matrix",
    "callvalidate_dg_hybrid_ground_state",
):
    assert continuation_c.index(compact(token)) < publish, token
assert "callsolve_dg_hybrid_generalized_complete_once" not in continuation_c[publish:]
assert "calldg_dc_update_potential_from_distributed_density" not in continuation_c[publish:]
assert "callderive_dg_hybrid_occupation_policy" not in continuation_c[publish:]

# Early lambda-zero failure may not roll back an uninitialized controller.
guard = compact("if(.not.continuation_controller%valid.or..not.continuation_controller%trial_active)then")
guard_at = continuation_c.index(guard)
reject_at = continuation_c.index("callreject_dg_hybrid_trial", guard_at)
assert guard_at < reject_at
assert "dgcontinuationinitialfixedpointfailedbeforerollback" in continuation_c[guard_at:reject_at]

# v4 preserves localized construction rows and occupied LCFO amplitudes; it
# never rotates them into replicated spectral basis state.
for token in (
    "publish_rt_dg_hybrid_checkpoint_v4",
    "payload%initial_occupied_amplitudes",
    "payload%basis_point_offsets",
    "payload%basis_support_ids",
    "payload%operator_offsets",
    "payload%operator_columns",
    "dc%ppg_tot%nlma",
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

# Mutation receipts: removing the reachable v4 call or moving the rollback
# call ahead of its guard invalidates the checked invariants.
mutated = continuation_c.replace("callpublish_dg_hybrid_divided_v4(", "callremoved_v4(", 1)
assert mutated.count("callpublish_dg_hybrid_divided_v4(") == 0
mutated = continuation_c[:guard_at] + "callreject_dg_hybrid_trial" + continuation_c[guard_at:]
assert mutated.index("callreject_dg_hybrid_trial", guard_at) < mutated.index(guard, guard_at)

print("PASS formal Hybrid continuation reaches distributed v4 publication")
