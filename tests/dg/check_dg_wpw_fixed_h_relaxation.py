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
assert "source_values(point,isource)=support_w_values(source_w_position)" in seed_body, (
    "source ensemble must be materialized from the same common-domain occupied-W descriptor used by the W basis"
)
assert "source_w_ids" in seed_body, (
    "source-column order must be mapped explicitly to occupied-W stable IDs"
)
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
assert "apply_dg_wpw_periodic_gradient_mpi" not in source_body, (
    "production must not require complete periodic-cell coverage"
)
assert "build_sawf_projected_buffer_gradients(occupied_stencil" in source_body, (
    "production must differentiate the projected fragment-box field"
)
assert "all(dc%nxyz_buffer>stencil_radius)" in source_body, (
    "user buffer must exceed the stencil radius and retain a DG tail layer"
)
assert "descriptor_buffer_lo=1-dc%nxyz_buffer+stencil_radius" in source_body and (
    "descriptor_buffer_hi=dc%nxyz_domain+dc%nxyz_buffer-stencil_radius" in source_body
), "the DG descriptor must be exactly the stencil-closed prepared-buffer domain"
assert "descriptor_gradient_lo=descriptor_buffer_lo" in source_body and (
    "descriptor_gradient_hi=descriptor_buffer_hi" in source_body
), (
    "DG values and gradients must use one common point set"
)
assert "gather_dg_wpw_occupied_w_payload" in source_body, (
    "spatially distributed occupied-W payload must be assembled by point ID"
)
gather_pos = source_body.index("call gather_dg_wpw_occupied_w_payload")
gather_call = source_body[gather_pos:gather_pos + 700]
assert "source_values(:,source_offset+1:source_offset+source_count)" in gather_call, (
    "descriptor core values must use the same direct core field as the density source"
)
assert ",source_core," not in gather_call, (
    "descriptor core values must retain stable-ID column ordering"
)
assert "reorder_dg_wpw_fragment_buffer" in source_body, (
    "fragment cyclic storage must be reordered into unwrapped P before differentiation"
)
assert "comm_summation(occupied_stencil_partial,occupied_storage" in source_body and (
    "info%icomm_r" in source_body
), "distributed storage pieces must be assembled before unwrapped P differentiation"
assert "any(mg%is>1)" not in source_body and "any(mg%ie<stencil_extent)" not in source_body, (
    "occupied-W P construction must not require every spatial rank to own the full fragment"
)
assert "descriptor_core_d_invariant" in source_body and "source_core(point,:)" in source_body, (
    "committed core values must be checked against direct D values after distributed reductions"
)
assert "source_values(core_point,source_index)=source_values(core_point,source_index)+received_value(record)" not in source_body, (
    "canonical aliases must not be added back into the already-periodic fragment field"
)
assert "initialize_dg_wpw_occupied_w_basis_collective" in source_body, (
    "projected-Wannier builder must transactionally commit the occupied-W descriptor"
)
for token in (
    "diagnose_sawf_discrete_wannier_spread(spread_link_global",
    "[dg-wpw-wannier-spread]",
    "omega_a2=",
    "width_a=",
    "center_valid=",
    "[dg-wpw-wannier-spread-summary]",
    "median_a=",
    "p90_a=",
    "count_above_1p2a=",
):
    assert token in source_body, f"missing occupied-W localization diagnostic: {token}"
spread_diagnostic_pos = source_body.index("call diagnose_sawf_discrete_wannier_spread")
buffer_gate_pos = source_body.index("'buffer_sufficiency'")
assert spread_diagnostic_pos < buffer_gate_pos, "spread diagnostic must precede the unchanged buffer gate"
assert "localize_sawf" not in source_body, "diagnostic phase must not rotate occupied W"
bootstrap_pos = LCFO.index("call assemble_wpw_core_p_bootstrap")
legacy_layout_pos = LCFO.index("call build_wpw_w_row_layout")
assert bootstrap_pos < legacy_layout_pos, (
    "W-independent core/P bootstrap must complete before legacy W layout during migration"
)
descriptor_build_pos = LCFO.index("call build_core_owned_projected_wannier_density_seed")
assert descriptor_build_pos < legacy_layout_pos, (
    "occupied-W descriptor must be committed before any W-row layout is created"
)
assert "build_dg_wpw_w_row_layout_from_owned_ids" in LCFO, (
    "WPW production W-row layout must use occupied-W descriptor IDs"
)
context_call = LCFO[LCFO.index("call initialize_dg_wpw_fragment_root_context"):]
assert "wpw_occupied_w_basis%global_count" in context_call[:500], (
    "production context W dimension must be the occupied-W global count, not legacy n_mat"
)
halo_start = LCFO.index("subroutine prepare_wpw_volume_halo")
halo_end = LCFO.index("end subroutine prepare_wpw_volume_halo", halo_start)
halo_body = LCFO[halo_start:halo_end]
assert "wpw_occupied_w_basis%buffer_values" in halo_body, (
    "WPW halo packing must read occupied-W descriptor buffer values"
)
assert "local_basis_value" not in halo_body and "basis_transform" not in halo_body, (
    "WPW halo packing must not reconstruct legacy Flux-LCFO rows"
)
quad_start = LCFO.index("subroutine assemble_wpw_core_volume_quadrature")
quad_end = LCFO.index("end subroutine assemble_wpw_core_volume_quadrature", quad_start)
quad_body = LCFO[quad_start:quad_end]
assert "evaluate_dg_wpw_occupied_w_point" in quad_body, (
    "WPW core quadrature must directly read the owned unwrapped D point"
)
assert "local_basis_value" not in quad_body and "local_basis_grad" not in quad_body, (
    "WPW core quadrature must not evaluate legacy Flux-LCFO W rows"
)
assert "if(info%id_o==0)then" in quad_body, (
    "orbital-replicated descriptor fields must be integrated on one orbital lane only"
)
assert "evaluate_dg_wpw_core_w_support" in quad_body and "add_dg_wpw_metric_point" in quad_body, (
    "core quadrature must integrate the halo-complete support-W Gram"
)
face_provider_start = LCFO.index("subroutine prepare_wpw_canonical_face_trace_provider")
face_provider_end = LCFO.index("end subroutine prepare_wpw_canonical_face_trace_provider", face_provider_start)
face_provider_body = LCFO[face_provider_start:face_provider_end]
assert "evaluate_dg_wpw_occupied_w_point" in face_provider_body, (
    "face traces must directly read the owned unwrapped D point"
)
assert "local_basis_value(" not in face_provider_body and "local_basis_grad(" not in face_provider_body, (
    "occupied-W face traces must not fall back to legacy positional LCFO rows"
)
assert "build_dg_wpw_metric_gram" in quad_body, (
    "fragment-rank WW/WP/PP Gram partials must reduce transactionally"
)
publish_start = LCFO.index("subroutine publish_wpw_core_volume_candidates")
publish_end = LCFO.index("end subroutine publish_wpw_core_volume_candidates", publish_start)
publish_body = LCFO[publish_start:publish_end]
assert "wpw_candidate_kind_ww" in publish_body and "target_fragment" in publish_body, (
    "real-space WW rows must route with explicit canonical owners"
)
assert "publish_dg_wpw_realspace_metric" in publish_body, (
    "routed WP/PP/WW Gram must be atomically published in the production context"
)
assert "owned%ww_kinetic" in publish_body and "owned_conjugate%ww_kinetic" in publish_body, (
    "routed WW kinetic entries and their conjugate audit must be atomically published"
)
assert "owned%ww_potential" in publish_body and "owned_conjugate%ww_potential" in publish_body, (
    "routed WW potential entries and their conjugate audit must be atomically published"
)
assert "wpw_metric_wp_h" in publish_body and "find_integer_id(wpw_owned_w_ids" not in publish_body, (
    "WP volume Hamiltonian must retain neighboring support-W tail contributions"
)
assert "saved_context" not in publish_body, "publication must not rely on allocating deep-copy rollback"
iter_start = LCFO.index("subroutine publish_wpw_iterated_operator")
iter_end = LCFO.index("end subroutine publish_wpw_iterated_operator", iter_start)
iter_body = LCFO[iter_start:iter_end]
assert "wpw_candidate_kind_ww" in iter_body and "root_ww_potential" in iter_body, (
    "dynamic WW potential must route complex support-W rows by stable ID"
)
assert "import_dg_wpw_lcfo_ww_components" not in iter_body, (
    "dynamic WW potential refresh must not rebuild positional LCFO components"
)
assert "ww_potential=old_ww_potential" in iter_body, (
    "dynamic WW potential rollback must restore WP/PP/WW atomically"
)
snapshot_start = LCFO.index("subroutine snapshot_wpw_fixed_point_state")
snapshot_end = LCFO.index("end subroutine restore_wpw_fixed_point_state", snapshot_start)
snapshot_body = LCFO[snapshot_start:snapshot_end]
assert "wpw_saved_ww_potential" in snapshot_body, (
    "fixed-point rollback must snapshot stable-ID real-space WW potential rows"
)
assert "wpw_saved_ww_components" not in snapshot_body and "publish_dg_wpw_lcfo_ww_components" not in snapshot_body, (
    "fixed-point rollback must not restore positional LCFO WW components"
)
face_start = LCFO.index("subroutine scan_wpw_canonical_faces")
face_end = LCFO.index("end subroutine scan_wpw_canonical_faces", face_start)
face_body = LCFO[face_start:face_end]
assert "support_w_owner" in face_body and "scan_dg_wpw_canonical_faces" in face_body, (
    "canonical WW face candidates must carry stable occupied-W row ownership"
)
projector_start = LCFO.index("subroutine assemble_wpw_projector_nonlocal")
projector_end = LCFO.index("end subroutine assemble_wpw_projector_nonlocal", projector_start)
projector_body = LCFO[projector_start:projector_end]
assert "wpw_metric_accumulator%w_points" in projector_body, (
    "projector overlaps must consume halo-complete descriptor support-W values"
)
assert "wpw_occupied_w_basis%global_count" in projector_body and "n_mat(1)" not in projector_body, (
    "projector W/P namespaces must use occupied-W provenance, not legacy LCFO row count"
)
assert "canonicalize_dg_wpw_projector_records" in projector_body and (
    "exchange_dg_wpw_projector_overlaps" not in projector_body
), "projector partials must reduce through canonical channel owners, not legacy fragment peer exchange"
assert "ppg%lma_tbl" in projector_body and "pp%nproj" in projector_body, (
    "stable projector keys must derive angular/projector identity from pseudopotential provenance"
)
assert "dg wpw occupied-w task 5 hamiltonian assembly is not yet installed" not in LCFO, (
    "completed occupied-W Hamiltonian route remains unconditionally stopped"
)
assert "call import_wpw_lcfo_ww_components" not in LCFO, (
    "production occupied-W route still mixes positional LCFO WW components"
)
for token in ("canonicalize_sawf_wannier_center", "canonicalize_sawf_bond_identity"):
    assert token in source_body, f"core-owned source builder does not connect {token}"
assert "exchange_sawf_discovered_wannier_tails" not in source_body, (
    "already-periodic P must not be canonicalized and accumulated a second time"
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
assert "wpw_bootstrap_source_values" in seed_body and "wpw_bootstrap_source_condition" in seed_body, (
    "density seed must reuse the pre-layout occupied-W bootstrap without rebuilding it"
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
