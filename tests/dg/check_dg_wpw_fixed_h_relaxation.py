#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GLOBAL = (ROOT / "src/io/salmon_global.f90").read_text().lower()
INPUT = (ROOT / "src/io/inputoutput.f90").read_text().lower()
LCFO = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text().lower()

control = "yn_dg_wpw_fixed_h_relaxation"
assert control in GLOBAL, "missing explicit fixed-H production control"
assert control in INPUT, "fixed-H control is absent from input/default/broadcast/log validation"
assert f"{control} = 'n'" in INPUT, "fixed-H mode must not change repository defaults"

for token in (
    "snapshot_wpw_frozen_h_state",
    "validate_wpw_frozen_h_state",
    "wpw_frozen_rho_tot",
    "wpw_frozen_vh_tot",
    "wpw_frozen_vxc_tot",
    "wpw_frozen_vloc_tot",
    "wpw_frozen_ww_components",
    "wpw_frozen_wp_volume",
    "wpw_frozen_pp_volume",
    "wpw_frozen_wp_nonlocal",
    "wpw_frozen_wp_face",
    "wpw_frozen_pp_nonlocal",
    "wpw_frozen_owned_w_ids",
    "wpw_frozen_required_w_ids",
    "wpw_frozen_owned_p_ids",
    "wpw_frozen_required_p_ids",
    "release_wpw_frozen_h_state",
    "wpw_frozen_state_local_matches",
):
    assert token in LCFO, f"missing frozen-H state contract: {token}"

assert "run_wpw_fixed_h_relaxation" in LCFO, "production route lacks fixed-H relaxation driver"
fixed_start = LCFO.index("subroutine run_wpw_fixed_h_relaxation")
fixed_end = LCFO.index("end subroutine run_wpw_fixed_h_relaxation", fixed_start)
fixed_body = LCFO[fixed_start:fixed_end]
assert "wpw_potential_step" not in fixed_body, "fixed-H relaxation must not feed reconstructed density into H"
assert "validate_wpw_frozen_h_state" in fixed_body, "fixed-H driver does not enforce frozen-state invariants"

snapshot_start = LCFO.index("subroutine snapshot_wpw_frozen_h_state")
snapshot_end = LCFO.index("end subroutine snapshot_wpw_frozen_h_state", snapshot_start)
snapshot_body = LCFO[snapshot_start:snapshot_end]
assert "wpw_frozen_ww_components=wpw_ww_components" not in snapshot_body, (
    "frozen snapshot must not rely on derived-type deep assignment"
)
assert "deallocate(wpw_frozen_ww_h0_dense," not in snapshot_body, (
    "frozen snapshot must not use grouped deallocation"
)
assert "stat=astat" in snapshot_body, "frozen snapshot allocations need stat= handling"

validate_start = LCFO.index("subroutine validate_wpw_frozen_h_state")
validate_end = LCFO.index("end subroutine validate_wpw_frozen_h_state", validate_start)
validate_body = LCFO[validate_start:validate_end]
assert "wpw_frozen_state_local_matches" in validate_body, (
    "frozen validation must be shape-safe before value comparisons"
)
assert "wpw_potential_stage_ok" in validate_body, (
    "frozen validation failures must be synchronized collectively"
)

seed_start = LCFO.index("subroutine build_wpw_projected_occupied_seed")
seed_end = LCFO.index("end subroutine build_wpw_projected_occupied_seed", seed_start)
seed_body = LCFO[seed_start:seed_end]
assert "root_overlap_p" in seed_body, "density-carrying seed must accumulate a P overlap block"
assert "wpw_volume_accumulator%p_points" in seed_body, "P overlap must use production P basis values"
assert "initialize_dg_wpw_metric_projected_occupied" in seed_body, (
    "density-carrying overlaps must be projected by solving S C=B"
)
assert "wpw_qp(:,1:wpw_nocc)=0" not in seed_body, "occupied P coefficients must not be zeroed"
assert "density_carrying_fragment_seed" in seed_body, "seed provenance must be explicit in diagnostics"
assert "coef_wf" not in seed_body, "density-carrying seed must not use Flux-LCFO eigenvectors"

print("PASS fixed-H WPW relaxation source contract")
