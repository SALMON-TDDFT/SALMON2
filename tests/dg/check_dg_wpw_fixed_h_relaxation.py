#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GLOBAL = (ROOT / "src/io/salmon_global.f90").read_text().lower()
INPUT = (ROOT / "src/io/inputoutput.f90").read_text().lower()
LCFO = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text().lower()
MATRIX_FREE = (ROOT / "src/gs/dc/dg_wpw_matrix_free_scf.f90").read_text().lower()

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
    "wpw_frozen_ww_projector_nonlocal",
    "wpw_frozen_ww_projector_cross_value",
    "wpw_frozen_ww_projector_cross_row_id",
    "wpw_frozen_ww_projector_cross_col_id",
    "wpw_frozen_callbacks_bound",
    "wpw_frozen_operator_epoch",
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
assert "set_dg_wpw_interface_lambda" in fixed_body, "fixed-H relaxation must start at zero interface"
assert "run_dg_wpw_matrix_free_algebra_step" in fixed_body, "fixed-H relaxation lacks a matrix-free eigensolve"
assert fixed_body.count("validate_wpw_frozen_h_state") >= 3, (
    "frozen state must be validated before seed/solve and after each solver stage"
)
assert "fixed_h_max_iter" in fixed_body, "fixed-H relaxation needs a bounded iteration cap"
assert "not yet installed" not in fixed_body, "fixed-H relaxation still contains the fail-closed placeholder"
assert "continue_wpw_fixed_h_interface" in fixed_body, "zero-interface solution is not continued to lambda one"

continuation_start = LCFO.index("subroutine continue_wpw_fixed_h_interface")
continuation_end = LCFO.index("end subroutine continue_wpw_fixed_h_interface", continuation_start)
continuation_body = LCFO[continuation_start:continuation_end]
for token in ("accepted_lambda", "trial_lambda", "lambda_step", "accepted_qw", "accepted_qp", "accepted_merit"):
    assert token in continuation_body, f"continuation transaction is missing {token}"
assert "lambda_step=0.5d0*lambda_step" in continuation_body, "rejected continuation steps must reduce step size"
assert "set_dg_wpw_interface_lambda" in continuation_body, "continuation cannot restore accepted lambda"
assert continuation_body.count("validate_wpw_frozen_h_state") >= 2, (
    "continuation must validate frozen state around every trial"
)

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

seed_start = LCFO.index("subroutine build_wpw_density_carrying_fragment_seed")
seed_end = LCFO.index("end subroutine build_wpw_density_carrying_fragment_seed", seed_start)
seed_body = LCFO[seed_start:seed_end]
source_start = LCFO.index("subroutine build_core_owned_projected_wannier_density_seed")
source_body = LCFO[source_start:seed_end]
source_end = LCFO.index("end subroutine build_core_owned_projected_wannier_density_seed", source_start)
source_builder_body = LCFO[source_start:source_end]
assert "wpw_core_p_accumulator" in source_builder_body, (
    "projected-Wannier builder must consume the W-independent core/P bootstrap"
)
assert "wpw_volume_accumulator" not in source_builder_body, (
    "projected-Wannier construction must not depend on legacy W storage"
)
assert "root_overlap_p" in seed_body, "density-carrying seed must accumulate a P overlap block"
assert "wpw_volume_accumulator%p_points" in seed_body, "P overlap must use production P basis values"
assert "evaluate_dg_wpw_core_w_support" in seed_body, "W overlap must include support rows on every fragment core"
assert seed_body.count("reconstruct_dg_wpw_core_w_support") >= 2, (
    "raw and normalized physical W values must retain neighboring support-W tails"
)
assert "[dg-wpw-metric-realspace-gram]" in seed_body, (
    "density-carrying seed must compare assembled S with the real-space Gram"
)
metric_gram_pos = seed_body.index("[dg-wpw-metric-realspace-gram]")
assert "call wpw_apply_s" in seed_body[:metric_gram_pos], (
    "metric/real-space diagnostic must explicitly apply production S to C_raw"
)
assert "[dg-wpw-assembled-metric-split]" in seed_body, (
    "assembled metric must be split into WW/WP/PP on C_raw"
)
assert "[dg-wpw-w-realspace-norm]" in seed_body, (
    "actual tail-carrying W functions need a global real-space norm probe"
)
assert "[dg-wpw-w-norm-local-max]" in seed_body, (
    "abnormal W norms must be split into owner-core and halo-tail contributions"
)
assert "[dg-wpw-w-halo-pack-max]" in LCFO, (
    "maximum abnormal W tail must be traced through halo packing"
)
assert "density_rw(iw,:)=density_rw_support(point_info,:)" not in seed_body, (
    "physical reconstruction must not discard nonowned support-W coefficients"
)
assert "reduce_dg_wpw_metric_rhs_partials" in seed_body, "W/P overlap partials must route to canonical owners"
assert "solve_dg_wpw_metric_projection" in seed_body, "density-carrying overlaps must solve S C_raw=B"
solve_pos = seed_body.index("call solve_dg_wpw_metric_projection")
capture_pos = seed_body.index("wpw_projection_captured_norm=")
normalize_pos = seed_body.index("call initialize_dg_wpw_projected_occupied", solve_pos)
assert solve_pos < capture_pos < normalize_pos, "captured norm must use C_raw before S-orthonormalization"
assert "wpw_occupations(iw)*real(capture(iw,iw),8)" in seed_body, (
    "captured norm must be occupation weighted"
)
assert "sum(wpw_occupations*source_norm_global)" in seed_body, (
    "captured norm must be relative to the occupation-weighted source norm"
)
assert "wpw_qp(:,1:wpw_nocc)=0" not in seed_body, "occupied P coefficients must not be zeroed"
assert "density_carrying_fragment_seed" in seed_body, "seed provenance must be explicit in diagnostics"
assert "coef_wf" not in seed_body, "density-carrying seed must not use Flux-LCFO eigenvectors"
assert "local_nocc=count(system%rocc" not in seed_body, (
    "core-owned source count must not count every occupied buffered-fragment eigenvector"
)
assert "if(info%id_o/=0)cycle" in seed_body, (
    "orbital-replicated projected Wannier values must contribute on one orbital lane only"
)
assert "build_sawf_projected_wannier_from_overlap" in source_body, (
    "core-owned source builder must reconstruct projected Wannier values from reduced overlaps"
)
assert "apply_sawf_projected_wannier_transform(occupied_buffer" in source_body, (
    "buffer values must reuse the core-projection polar transform and Wannier gauge"
)
for token in ("canonicalize_sawf_wannier_center", "canonicalize_sawf_bond_identity"):
    assert token in source_body, f"core-owned source builder does not connect {token}"
assert "exchange_sawf_discovered_wannier_tails" in source_body, (
    "core-owned source builder must exchange dynamically discovered buffer tails"
)
assert "qualify_sawf_wannier_buffer_tail" in source_body, (
    "core-owned source builder must reject an insufficient outer buffer shell"
)
assert source_body.count("wpw_seed_collective_stage_ok") >= 5, (
    "Wannier source phases must synchronize local failures before later collectives"
)
for token in (
    "build_sawf_wannier_density",
    "transform_sawf_wannier_occupation",
    "qualify_sawf_wannier_density_projection",
):
    assert token in seed_body, f"production seed does not connect {token}"
assert "core_owned_projected_wannier_density_seed" in seed_body, (
    "checkpoint diagnostics lack projected-Wannier source provenance"
)
assert "optional,intent(in)::stagnation_limit" in MATRIX_FREE, (
    "metric solver does not expose a bounded diagnostic stagnation window"
)
assert "stagnation_limit=257" in seed_body, (
    "Si64 diagnostic route does not distinguish early stagnation from max-iteration nonconvergence"
)
assert "optional,intent(in)::diagnose_recurrence" in MATRIX_FREE, (
    "metric solver cannot compare recursive and explicitly recomputed residuals"
)
assert "diagnose_recurrence=.true." in seed_body, (
    "Si64 metric diagnostic does not enable recurrence/Hermitian checks"
)
for token in ("dg-wpw-metric-ritz", "near_null_weight", "ritz_relative_min"):
    assert token in MATRIX_FREE, f"metric solver lacks PCG/Lanczos spectral diagnostic: {token}"

print("PASS fixed-H WPW relaxation source contract")
