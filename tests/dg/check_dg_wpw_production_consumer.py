#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
MAIN = ROOT / "src/gs/main_dft.f90"
LCFO = ROOT / "src/gs/dc/lcfo_flux.f90"
CONTEXT = ROOT / "src/common/dg_wpw_production_context.f90"
BUILDER = ROOT / "src/rt/dg/rt_dg_wpw_production_builder.f90"
SCANNER = ROOT / "src/rt/dg/rt_dg_wpw_face_trace_scanner.f90"
BOUNDED = ROOT / "src/common/dg_wpw_bounded_operator.f90"
COMPLEMENT = ROOT / "src/common/dg_wpw_s_orthogonal_complement.f90"

main = MAIN.read_text().lower()
lcfo = LCFO.read_text().lower()
assert CONTEXT.exists(), "missing GS-neutral production context"
assert BUILDER.exists(), "missing rank-local production operator builder"

context = CONTEXT.read_text().lower()
builder = BUILDER.read_text().lower()
scanner = SCANNER.read_text().lower()
owner_exchange = (ROOT / "src/common/dg_wpw_owner_exchange.f90").read_text().lower()
bounded = BOUNDED.read_text().lower()
assert COMPLEMENT.exists(), "missing distributed S-orthogonal PW complement"
complement = COMPLEMENT.read_text().lower()

for token in (
    "type,public::s_dg_wpw_s_orthogonal_complement",
    "initialize_dg_wpw_s_orthogonal_complement",
    "apply_h_dg_wpw_s_orthogonal_complement",
    "apply_s_dg_wpw_s_orthogonal_complement",
    "map_dg_wpw_complement_to_original",
    "map_dg_wpw_original_to_complement",
    "validate_dg_wpw_s_orthogonal_complement",
    "release_dg_wpw_s_orthogonal_complement",
):
    assert token in complement, f"missing complement API: {token}"
for forbidden in ("global_sww", "global_a", "root_sww", "root_a"):
    assert forbidden not in complement, f"single-rank global complement allocation: {forbidden}"
for required in ("fetch_rows_from_owners", "reduce_w_partial_to_owners", "global_gram_dg_wpw_bounded"):
    assert required in complement, f"distributed complement must use {required}"
for required in ("solve_sww_deflated_block", "solve_sww_block_pcg", "diagnose_sww", "diagnose_sperp"):
    assert required in complement, f"missing distributed complement algorithm: {required}"
assert "m=size(global_w_ids)" in complement.replace(" ", ""), \
    "S_WW cutoff certification must span the full occupied-W dimension"
assert complement.count("call mpi_allreduce(local_bad,global_bad") >= 6, \
    "distributed solve/actions require collective validation"

adapter = "apply_dg_wpw_fragment_block_eigen_correction"
assert f"public::{adapter}" in bounded, "metric block eigensolver adapter is not public"
adapter_start = bounded.index(f"subroutine {adapter}")
adapter_end = bounded.index("end subroutine", adapter_start)
adapter_body = bounded[adapter_start:adapter_end]
collective = adapter_body.index("call mpi_allreduce(local_bad,global_bad")
delegate = adapter_body.index("call apply_dg_wpw_fragment_block_preconditioner")
assert collective < delegate, "adapter must validate metadata collectively before delegation"
assert "ieee_is_finite(eigenvalues)" in adapter_body and "size(eigenvalues)/=size(rw,2)" in adapter_body
assert "op%w_schedule%comm" in adapter_body

for token in (
    "initialize_dg_wpw_fragment_root_context",
    "owned_w_ids",
    "support_w_ids",
    "build_dg_wpw_rank_local_quadrature",
    "scan_dg_wpw_canonical_faces",
    "build_windowed_sparse_wpw_operators",
    "bind_dg_wpw_hs_callbacks",
    "consume_dg_wpw_bounded_subspace",
):
    assert token in context + builder + main, f"missing production integration point: {token}"

for text, label in ((context, "context"), (builder, "builder")):
    for name in ("h_global", "s_global", "wp_global", "pp_global"):
        assert f"allocate({name}" not in text, \
            f"global dense operator allocation in production {label}: {name}"

assert not re.search(r"do\s+ifrag\s*=\s*1\s*,\s*n_frag", scanner), \
    "canonical scanner must not traverse all fragments"
assert not re.search(r"\bmpi_(?:allreduce|allgather|send|recv|isend|irecv)\b", scanner), \
    "canonical scanner contains hidden MPI"
assert "pp_face" not in builder, "periodic-H1 PP must not have face candidates"
assert "(k-1)*n_g+g_id" in context.replace(" ", ""), \
    "stable windowed_kg column id contract is missing"

for token in (
    "select_dg_wpw_g_modes",
    "initialize_dg_wpw_fragment_root_context",
    "build_dg_wpw_rank_local_quadrature",
    "route_dg_wpw_staged_candidates",
    "prepare_dg_wpw_trace_halo",
    "scan_dg_wpw_canonical_faces",
    "build_dg_wpw_production_operator",
    "bind_dg_wpw_hs_callbacks",
    "run_wpw_production_scf",
    "advance_dg_wpw_scf_iteration",
    "verify_dg_wpw_scf_fixed_point",
    "fetch_dg_wpw_support_coefficients",
    "build_dg_wpw_core_density",
    "update_wpw_owned_lda_hartree",
    "evaluate_and_run_dg_wpw_lcfo_potential_map",
    "route_and_replace_dg_wpw_potential_volume",
):
    assert token in lcfo, f"lcfo_flux does not publish production WPW state: {token}"

assert re.search(r"call\s+run_wpw_production_scf\b", lcfo), \
    "lcfo_flux does not execute the stepwise production SCF after callback binding"
assert "wpw_fixed_point_mode" in lcfo and re.search(
    r"merge\s*\(\s*1d0\s*,\s*dg_wpw_scf_mix\s*,\s*wpw_fixed_point_mode\s*\)", lcfo
), "production fixed-point verification must rebuild with density mixing fixed to one"
assert "dc%system_tot%rocc" not in lcfo, \
    "DCDft deallocates total-system rocc before LCFO; production must not dereference it"
assert "remaining_electrons=dc%elec_num_tot" in lcfo.replace(" ", "") and \
       "min(2d0,remaining_electrons)" in lcfo.replace(" ", ""), \
    "production must construct sharp nspin=1 occupations from the total electron count"
assert re.search(r"mpi_bcast\s*\(\s*wpw_eigenvalues\b", lcfo), \
    "fragment nonroots must receive production eigenvalues before evaluating total energy"
assert "potential_norms_local" in lcfo and re.search(
    r"mpi_allreduce\s*\(\s*potential_norms_local\s*,\s*potential_norms_global", lcfo
), "potential residual must be globally reduced before the SCF driver compares scalar values"
assert lcfo.count("wpw_potential_stage_ok(stage_info)") >= 4, \
    "rank-local potential failures must be reduced before peers enter later collectives"
assert "snapshot_wpw_fixed_point_state" in lcfo and "restore_wpw_fixed_point_state" in lcfo, \
    "failed fixed-point probes must restore operator, accumulator, and total-grid fields"
assert "energy_min_total" in lcfo and "energy_max_total" in lcfo, \
    "candidate energy must be checked across all fragments, not only within each fragment communicator"
assert "old_wp_volume" in lcfo and "old_pp_volume" in lcfo and \
       lcfo.count("replace_dg_wpw_potential_volume(wpw_context,old_wp_volume,old_pp_volume") >= 3, \
    "failed potential publication/rebuild/binding must restore the prior valid WP/PP volume"
assert "if(.not.sawf_explicit_basis_active.or..not.allocated(sawf_explicit_buffer))" not in \
       lcfo.replace(" ", ""), \
    "production volume halos must allow the LCFO basis_transform+spsi buffered fallback"
assert "wpw_window_buffer=merge(dc%nxyz_buffer-stencil_radius,0,num_fragment>1)" in lcfo.replace(" ", ""), \
    "production windows must use transition buffers only along fragmented axes"
assert "wpw_trace_face_coord" in lcfo and "dc%nxyz_domain(axis)+1" in lcfo.replace(" ", ""), \
    "plus face trace must use the first buffer point shared with the neighbor core"

volume = lcfo[lcfo.index("subroutine assemble_wpw_core_volume_quadrature"):
              lcfo.index("end subroutine assemble_wpw_core_volume_quadrature")]
assert "size(wpw_owned_w_ids)" in volume and \
       "add_dg_wpw_core_point(wpw_volume_accumulator,owned_w,owned_grad_w" in volume.replace(" ", ""), \
    "piecewise DG volume assembly must use fragment-owned W rows only"
assert "evaluate_dg_wpw_core_w_support" in volume and \
       "add_dg_wpw_metric_point(wpw_metric_accumulator" in volume.replace(" ", ""), \
    "tail-carrying support W must feed the real-space metric accumulator"

potential = lcfo[lcfo.index("subroutine wpw_potential_step"):
                 lcfo.index("end subroutine wpw_potential_step")]
assert "rw(size(wpw_owned_w_ids),wpw_nocc)" in potential.replace(" ", ""), \
    "SCF density must use owned-W coefficients consistently with orthonormal_ww"

faces = lcfo[lcfo.index("subroutine scan_wpw_canonical_faces"):
             lcfo.index("end subroutine scan_wpw_canonical_faces")]
assert "wpw_support_w_ids" in faces, \
    "canonical face assembly must retain neighbor W support"

assert re.search(r"\bcall\s+build_dg_wpw_rank_local_quadrature\b", lcfo), \
    "lcfo_flux does not execute fragment-rank WPW volume quadrature"
assert "projector_weights=ppg%rinv_uvu(1:nlma)/hvol" in lcfo.replace(" ", ""), \
    "weighted projector overlaps require rinv_uvu/hvol to match SALMON dpseudo's single hvol convention"
for stage in ("canonical_projector_owner", "ww_local", "ww_cross", "wp", "pp", "install"):
    assert f"call trace_wpw_projector_stage('{stage}'" in lcfo, \
        f"missing bounded projector-stage diagnostic: {stage}"
assert "size(support_records)" in lcfo and "size(wpw_context%pp_r)" in lcfo, \
    "projector diagnostics must expose record and PP cardinalities before costly assembly"
for stage in ("candidate_copy", "sparse_blocks", "bounded_operator", "block_publication"):
    assert f"call trace_production_builder('{stage}'" in builder, \
        f"missing production-builder boundary diagnostic: {stage}"
for stage in ("graph", "id_exchange", "owner_resolution"):
    assert f"call trace_owner_schedule('{stage}'" in owner_exchange, \
        f"missing owner-schedule boundary diagnostic: {stage}"
for stage in ("dimensions", "seed", "state", "algebra_begin", "algebra_end"):
    assert f"call trace_wpw_scf_stage('{stage}'" in lcfo, \
        f"missing production-SCF boundary diagnostic: {stage}"

assert "yn_dg_wpw_production" in main, "production route is not explicit in main_dft"
assert "run_dg_wpw_fixed_scf" not in main, "dense oracle must be unreachable from main_dft"

print("PASS distributed WPW production consumer source contract")
