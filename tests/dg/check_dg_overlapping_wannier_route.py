#!/usr/bin/env python3
"""Source contract for the isolated overlapping-Wannier DC route."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

for retained_contract_path in (
    "tests/dg/check_obsolete_dg_routes_removed.py",
    "docs/plans/2026-07-31-obsolete-dg-route-inventory.md",
    "src/gs/dc/dg_overlapping_wannier_checkpoint.f90",
    "src/gs/dc/dg_overlapping_wannier_construction.f90",
    "src/gs/dc/dg_overlapping_wannier_localization.f90",
    "src/gs/dc/dg_overlapping_wannier_w90.f90",
    "src/rt/dg/rt_dg_overlapping_wannier.f90",
    "src/gs/dc/lcfo.f90",
    "src/gs/eigen_subdiag_eigenexa.f90",
):
    assert (ROOT / retained_contract_path).is_file(), (
        f"missing retained route contract/source: {retained_contract_path}"
    )


def source(path: str) -> str:
    return (ROOT / path).read_text()


global_source = source("src/io/salmon_global.f90")
input_source = source("src/io/inputoutput.f90")
main_source = source("src/gs/main_dft.f90")
eigenexa_source = source("src/gs/eigen_subdiag_eigenexa.f90")
rt_main_source = source("src/rt/main_tddft.f90")
scf_source = source("src/gs/scf_iteration_dft.f90")
dcdft_source = source("src/gs/dc/dcdft.f90")
types_source = source("src/gs/dc/dg_overlapping_wannier_types.f90")
construction_source = source("src/gs/dc/dg_overlapping_wannier_construction.f90")
assert construction_source.lower().count(
    "product_table(right_operation,left_operation)"
) >= 2, "point pullback composition must use the reversed geometric product order"
localization_source = source("src/gs/dc/dg_overlapping_wannier_localization.f90")
w90_source = source("src/gs/dc/dg_overlapping_wannier_w90.f90")
assert re.match(r"\s*#include\s+[\"<]config\.h[\">]", w90_source), (
    "Wannier90 adapter must import CMake feature macros before conditional compilation"
)
lcfo_source = source("src/gs/dc/lcfo.f90")
projection_source = source("src/gs/dc/dg_overlapping_wannier_projection.f90")
operators_source = source("src/gs/dc/dg_overlapping_wannier_operators.f90")
ow_scf_source = source("src/gs/dc/dg_overlapping_wannier_scf.f90")
ow_solver_source = source("src/gs/dc/dg_overlapping_wannier_solver.f90")
ow_checkpoint_source = source("src/gs/dc/dg_overlapping_wannier_checkpoint.f90")
ow_rt_source = source("src/rt/dg/rt_dg_overlapping_wannier.f90")
dc_cmake = source("src/gs/dc/CMakeLists.txt")
w90_builder_source = source("cmakefiles/Builder/build_wannier90.cmake")
xc_source = source("src/xc/salmon_xc.f90")
si64_runner_source = source("tests/dg/run_si64_overlapping_wannier_gate.py")
si64_checker_source = source("tests/dg/check_si64_overlapping_wannier_gate.py")

ow_ground_state = re.search(
    r"subroutine\s+run_dg_overlapping_wannier_ground_state_for_main(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert ow_ground_state
ow_ground_state_body = ow_ground_state.group("body").lower()
for mlwf_call in (
    "assemble_dg_w90_gamma_matrices",
    "run_dg_w90_gamma_library",
    "apply_dg_w90_gamma_transform",
):
    assert re.search(rf"call\s+{mlwf_call}\b", ow_ground_state_body), (
        f"production overlapping-Wannier route must call {mlwf_call}"
    )
assert "call localize_dg_occupation_blocks" not in ow_ground_state_body, (
    "production overlapping-Wannier V3 route must not call the custom localizer"
)
for provenance_field in (
    "mlwf_backend",
    "mlwf_version",
    "mlwf_input_fingerprint",
    "mlwf_transform_fingerprint",
    "mlwf_spreads",
    "mlwf_coordinator_bytes",
    "mlwf_workspace_peak_bytes",
    "mlwf_coordinator_byte_limit",
    "mlwf_symmetry_receipts",
    "mlwf_canonical",
    "affine_group_order",
    "translation_subgroup_order",
    "point_cogroup_order",
    "fixed_center_group_order",
    "fixed_center_group_fingerprint",
    "fixed_center_fractional",
    "fixed_center_inversion_present",
    "affine_proof_workspace_peak_bytes",
    "point_projection_workspace_peak_bytes",
    "occupied_subspace_distance",
    "occupied_electron_count_drift",
    "occupied_density_interior_difference",
    "occupied_density_boundary_difference",
    "occupied_density_interior_tolerance",
    "occupied_density_boundary_tolerance",
    "occupied_closure_before",
    "occupied_closure_after",
    "occupied_selected_edge",
    "occupied_rejected_edge",
    "occupied_cluster_gap",
    "occupied_selected_block_dimension",
    "occupied_adaptation_workspace_peak_bytes",
):
    assert provenance_field in ow_checkpoint_source.lower(), (
        f"V3 checkpoint must serialize and validate {provenance_field}"
    )
    assert re.search(
        rf"ow_checkpoint%{provenance_field}\s*=", main_source, re.I
    ), f"production publication must bind {provenance_field}"
assert "fingerprint_ow_w90_matrices" in ow_ground_state_body
assert "fingerprint_ow_w90_transform" in ow_ground_state_body
assert re.search(
    r"call\s+orthonormalize_dg_distributed_seed_space\b", ow_ground_state_body
), "affine-closed LCFO seeds must be orthonormalized without full-orbit regeneration"
adaptation_call = ow_ground_state_body.find(
    "call build_dg_group_averaged_occupied_candidates_eigenexa"
)
affine_measurement_call = ow_ground_state_body.rfind(
    "call measure_dg_rank_fixed_symmetry_residuals_eigenexa"
)
dmn_begin_call = ow_ground_state_body.find("call begin_sawf_dmn")
assert 0 <= adaptation_call < affine_measurement_call < dmn_begin_call, (
    "fixed-center group-averaged occupied selection must precede full-affine proof and DMN"
)
assert re.search(
    r"global_seed_values\s*\(\s*1\s*:\s*nstate\s*,\s*:\s*\)\s*=\s*"
    r"(?:adapted_)?occupied_candidates",
    ow_ground_state_body,
), "production must replace only the occupied seed block by symmetry-adapted candidates"
assert not re.search(
    r"call\s+build_dg_distributed_symmetry_closed_basis\b", ow_ground_state_body
), "production must not regenerate every affine image after its closure receipt passes"
for forbidden_dense_affine in (
    "global_symmetry_overlap",
    "global_candidate_raw",
    "global_candidate_representation",
    "global_retained_representation",
):
    assert forbidden_dense_affine not in ow_ground_state_body, (
        f"production must not retain full-affine dense tensor {forbidden_dense_affine}"
    )
assert "lcfo_symmetry_worst_operation" in ow_ground_state_body
assert not re.search(
    r"do\s+io\s*=\s*1\s*,\s*size\s*\(\s*global_point_product.*?"
    r"LCFO_symmetry_operation=",
    ow_ground_state_body,
    re.I | re.S,
), "production must summarize full-affine residuals instead of logging every operation"

assert re.search(r"call\s+zheev\s*\(\s*'v'\s*,\s*'u'", ow_solver_source, re.I), (
    "the bounded reduced Hermitian Ritz problem must use the existing LAPACK path"
)
assert "do sweep=1,100" not in ow_solver_source, (
    "production block Ritz must not use the replicated cubic Jacobi sweep loop"
)
assert "call append_s_orthonormal_block" in ow_solver_source, (
    "block residual expansion must batch its metric projection and reduction"
)
assert "0.2d0*density_tolerance" in ow_scf_source, (
    "outer SCF must leave margin for the mandatory unmixed fixed-point gate"
)
for text in (si64_runner_source, si64_checker_source):
    assert '"buffer5"' in text and '"buffer6"' in text, (
        "Si64 matrix must compare admissible split-axis buffer depths"
    )
assert "(5, 5, 5)" in si64_runner_source and "(6, 6, 6)" in si64_runner_source
assert "4x2x1" not in si64_runner_source and "4x2x1" not in si64_checker_source, (
    "the 32-cubed Si64 acceptance matrix must keep the physical decomposition fixed"
)
assert si64_checker_source.count("validate_fixed_decomposition(variables") >= 2, (
    "every Si64 acceptance row and the normal reference must verify runtime num_fragment"
)
assert 'PRODUCTION_BOX_PROFILE = "buffer5"' in si64_checker_source, (
    "Task 9 must name buffer5 explicitly as the fixed production profile"
)
assert "itertools.combinations(values, 2)" not in si64_checker_source, (
    "buffer6 is diagnostic evidence, not a cross-buffer acceptance coordinate"
)
assert re.search(
    r"args\.minimum_row.*?\(\s*\"2x2x2\"\s*,\s*\"buffer6\"\s*,\s*\"c192-sp48\"",
    si64_runner_source,
    re.S,
), "the focused smoke gate must use the demonstrated buffer-6 configuration"
assert "manifest_digest_before" in si64_runner_source
assert "route checkpoint changed during restart reuse" in si64_runner_source
assert "SALMON_OW_GS_CHECKPOINT_V3" in si64_checker_source
assert "SALMON_OW_GS_RANK_SHARD_V3" in si64_checker_source
for token in ("validate_restart_log", "[ow-scf-diagnostic]", "forbidden route marker"):
    assert token in si64_checker_source.lower(), (
        "restart validation must reject recomputation and forbidden route markers"
    )
assert r"=\s*(\S+)" in si64_checker_source, (
    "Si64 evidence parser must accept formatted Fortran whitespace after '='"
)
assert re.search(
    r"yn_dg_dc_overlapping_wannier\s*==\s*'y'.*?"
    r"any\s*\(\s*num_rgrid\s*/\s*num_fragment\s*\+\s*2\s*\*\s*num_rgrid_buffer\s*>\s*num_rgrid\s*\)",
    input_source,
    re.I | re.S,
), "overlapping-Wannier buffer box must not exceed the periodic system on any axis"
assert re.search(
    r"yn_dg_dc_overlapping_wannier\s*==\s*'y'.*?"
    r"all\s*\(\s*num_rgrid\s*/\s*num_fragment\s*\+\s*2\s*\*\s*num_rgrid_buffer\s*==\s*num_rgrid\s*\)",
    input_source,
    re.I | re.S,
), "overlapping-Wannier buffer box must not become the complete periodic system"
assert "maximum_block_size=min(n,3*nstate)" in ow_solver_source, (
    "thick-restart coefficient iteration must bound its X/R/P space by n"
)
assert "previous_direction=old_q" in ow_solver_source, (
    "thick restart must retain the preceding Ritz subspace"
)
assert re.search(
    r"block_size\s*==\s*n.*?coefficient_diagnostics.*?"
    r"full-space Ritz residual exceeds numerical quality gate",
    ow_solver_source,
    re.I | re.S,
), (
    "a full-space Ritz solve must use the final numerical-quality gate instead of "
    "constructing another rank-deficient search space"
)
append_body = ow_solver_source[
    ow_solver_source.index("subroutine append_s_orthonormal_block"):
    ow_solver_source.index("end subroutine", ow_solver_source.index("subroutine append_s_orthonormal_block"))
]
assert re.search(r"do\s+projection_pass\s*=\s*1\s*,\s*2", append_body, re.I), (
    "small residual directions need two-pass S projection to avoid cancellation-amplified "
    "loss of orthogonality"
)

flag = r"yn_dg_dc_overlapping_wannier"
rt_flag = r"yn_dg_overlapping_wannier_rt"
rt_restart_flag = r"yn_dg_overlapping_wannier_rt_restart"

assert re.search(rf"character\s*\(\s*1\s*\).*::\s*{rt_flag}", global_source, re.I)
assert re.search(rf"\b{rt_flag}\s*=\s*'n'", input_source, re.I)
assert re.search(rf"namelist\s*/\s*propagation\s*/.*?\b{rt_flag}\b", input_source, re.I | re.S)
assert re.search(rf"call\s+comm_bcast\s*\(\s*{rt_flag}\b", input_source, re.I)
assert re.search(rf"write\s*\(\s*fh_variables_log.*?{rt_flag}", input_source, re.I | re.S)
assert re.search(rf"call\s+yn_argument_check\s*\(\s*{rt_flag}\s*\)", input_source, re.I)
assert re.search(rf"character\s*\(\s*1\s*\).*::\s*{rt_restart_flag}", global_source, re.I)
assert re.search(rf"\b{rt_restart_flag}\s*=\s*'n'", input_source, re.I)
assert re.search(rf"call\s+yn_argument_check\s*\(\s*{rt_restart_flag}\s*\)", input_source, re.I)
assert re.search(r"if\s*\(\s*yn_restart\s*==\s*'y'\s*\).*?forbids conventional", input_source, re.I | re.S)
assert re.search(
    rf"if\s*\(\s*{rt_flag}\s*==\s*'y'\s*\)\s*then.*?"
    r"call\s+run_dg_overlapping_wannier_coefficient_rt.*?return",
    rt_main_source,
    re.I | re.S,
), "coefficient RT needs a terminating dispatch before conventional RT initialization"
rt_dispatch = rt_main_source.lower().find("call run_dg_overlapping_wannier_coefficient_rt")
legacy_init = rt_main_source.lower().find("call initialization_rt")
assert 0 <= rt_dispatch < legacy_init
coefficient_driver = re.search(
    r"subroutine\s+run_dg_overlapping_wannier_coefficient_rt(?P<body>.*?)end\s+subroutine",
    rt_main_source,
    re.I | re.S,
)
assert coefficient_driver
coefficient_body = coefficient_driver.group("body")
compact_coefficient_body = re.sub(r"\s+", "", coefficient_body.lower())
for token in (
    "step==1.and.state%step==0.and.trim(ae_shape1)=='impulse'",
    "electric_field=-vector_potential_samples(:,step)/dt",
):
    assert token in compact_coefficient_body, (
        f"coefficient RT must preserve the fresh impulse jump: {token}"
    )
for required in (
    "read_dg_overlapping_wannier_checkpoint",
    "initialize_dg_overlapping_wannier_rt",
    "advance_dg_overlapping_wannier_rt",
    "read_dg_overlapping_wannier_rt_restart",
    "write_dg_overlapping_wannier_rt_restart",
    "calc_ac_ext_t",
):
    assert required in coefficient_body.lower()
for forbidden in (
    "initialization_rt",
    "time_evolution_step",
    "time_evolution_dg_fragment",
    "dc_lcfo",
    "eigenexa",
    "dg_wpw",
    "checkpoint_rt",
):
    assert forbidden not in coefficient_body.lower(), (
        f"coefficient RT must not enter forbidden route: {forbidden}"
    )
assert "close(unit,iostat=close_ios)" in ow_rt_source.lower()
assert "stored_digest/=rt_restart_digest" in ow_rt_source.lower()
assert "call zhegv" in ow_rt_source.lower()
assert "crank" not in ow_rt_source.lower()
for observable_api in (
    "evaluate_dg_overlapping_wannier_observables",
    "write_dg_overlapping_wannier_rt_observable_sample",
):
    assert re.search(rf"public\s*::.*?{observable_api}", ow_rt_source, re.I | re.S), (
        f"coefficient RT must expose {observable_api}"
    )
for observable_call in (
    "evaluate_dg_overlapping_wannier_observables",
    "write_dg_overlapping_wannier_rt_observable_sample",
):
    assert re.search(rf"call\s+{observable_call}", coefficient_body, re.I), (
        f"production coefficient RT must call {observable_call}"
    )
assert "checkpoint%occupations" in coefficient_body.lower(), (
    "production observables must use V3 checkpoint occupations"
)
assert "overlapping_wannier_rt_observables.dat" in coefficient_body.lower(), (
    "production coefficient RT must publish the dedicated observable time series"
)
assert "step time ex ey ez px py pz jx jy jz" in ow_rt_source.lower(), (
    "observable evidence must declare the deterministic E/P/J column order"
)

assert re.search(rf"character\s*\(\s*1\s*\).*::\s*{flag}", global_source, re.I), (
    "the overlapping-Wannier route needs its own global flag"
)
assert re.search(rf"\b{flag}\s*=\s*'n'", input_source, re.I), (
    "the new route must default off"
)
assert re.search(rf"namelist\s*/\s*dc\s*/.*?\b{flag}\b", input_source, re.I | re.S), (
    "the new flag must be part of the DC namelist"
)
assert re.search(rf"call\s+comm_bcast\s*\(\s*{flag}\b", input_source, re.I), (
    "the new flag must be broadcast collectively"
)
assert re.search(
    rf"write\s*\(\s*fh_variables_log.*?{flag}", input_source, re.I | re.S
), "the selected route must be recorded in variables.log"
assert re.search(rf"call\s+yn_argument_check\s*\(\s*{flag}\s*\)", input_source, re.I), (
    "the new flag must accept only y/n"
)
for tolerance in (
    "dg_ow_boundary_value_tolerance",
    "dg_ow_boundary_gradient_tolerance",
    "dg_ow_symmetry_tolerance",
    "dg_ow_localization_support_tolerance",
    "dg_ow_localization_spread_tolerance",
    "dg_ow_localization_gradient_tolerance",
):
    assert tolerance in global_source
    assert tolerance in input_source
assert re.search(r"dg_ow_localization_support_tolerance\s*=\s*1d-3", input_source, re.I), (
    "default localization graph must remain sparse for production-sized retained blocks"
)
assert re.search(r"dg_ow_localization_spread_tolerance\s*=\s*1d-14", input_source, re.I), (
    "localization spread resolution must retain resolvable double-precision descent"
)
assert re.search(r"dg_ow_localization_max_iterations\s*=\s*1024", input_source, re.I), (
    "localization iteration budget must cover the genuine-Si64 RCG convergence envelope"
)
assert "dg_ow_localization_max_iterations" in global_source
assert "dg_ow_localization_max_iterations" in input_source
si64_gs_input = (ROOT / "tests/dg/data/si64_overlapping_wannier_rt/input_gs.in").read_text()
assert re.search(r"dg_ow_localization_gradient_tolerance\s*=\s*2d-2", si64_gs_input, re.I), (
    "qualitative Si64 evidence must use the reviewed 0.02 absolute spread-gradient gate"
)
for window in (
    "dg_ow_candidate_states_per_fragment",
    "dg_ow_target_wanniers_per_fragment",
):
    assert window in global_source
    assert window in input_source

route_checks = [
    (r"trim\s*\(\s*theory\s*\)\s*/=\s*'dft'", "ground-state DFT only"),
    (r"\byn_dc\s*/=\s*'y'", "DC only"),
    (r"\byn_periodic\s*/=\s*'y'", "periodic only"),
    (r"\byn_spinorbit\s*==\s*'y'", "non-SOI only"),
    (
        r"num_kgrid\s*\(\s*1\s*\)\s*\*\s*num_kgrid\s*\(\s*2\s*\)\s*\*\s*"
        r"num_kgrid\s*\(\s*3\s*\)\s*/=\s*1",
        "Gamma only",
    ),
    (r"trim\s*\(\s*xc\s*\)\s*/=\s*'pz'", "PZ LDA only"),
    (r"\byn_dc_lcfo\s*==\s*'y'", "LCFO forbidden"),
    (r"\byn_self_checkpoint\s*==\s*'y'", "normal checkpoint forbidden"),
    (r"\bcheckpoint_interval\s*>=\s*1", "periodic checkpoint forbidden"),
]
for condition, requirement in route_checks:
    assert re.search(
        rf"if\s*\(\s*{flag}\s*==\s*'y'.*?{condition}",
        input_source,
        re.I | re.S,
    ), f"overlapping-Wannier validation must enforce: {requirement}"

assert not re.search(
    rf"if\s*\(\s*{flag}\s*==\s*'y'.*?\byn_eigenexa\s*==\s*'y'.*?forbids EigenExa",
    input_source,
    re.I | re.S,
), "the internal LCFO coefficient path must permit EigenExa"
assert re.search(r"yn_eigenexa\s*=\s*'y'", si64_gs_input, re.I), (
    "strict Si64 must exercise the LCFO EigenExa coefficient path"
)

assert re.search(
    rf"if\s*\(\s*{flag}\s*==\s*'y'\s*\)\s*then.*?"
    r"call\s+run_dg_overlapping_wannier_ground_state_for_main.*?"
    r"(?:return|error\s+stop)",
    main_source,
    re.I | re.S,
), "main_dft needs an explicit, terminating new-route dispatch"
scf_position = main_source.lower().find("call scf_iteration_dft")
route_dispatch_position = main_source.lower().find(
    "call run_dg_overlapping_wannier_ground_state_for_main"
)
lcfo_position = main_source.lower().find("call dc_lcfo_flux", route_dispatch_position)
assert 0 <= scf_position < route_dispatch_position < lcfo_position, (
    "construction dispatch must consume the conventional candidate window "
    "after SCF and terminate before LCFO"
)
assert "register_dg_overlapping_wannier_route_driver" not in ow_scf_source
assert "execute_registered_dg_overlapping_wannier_ground_state" not in ow_scf_source
assert re.search(
    r"subroutine\s+run_dg_overlapping_wannier_ground_state_for_main.*?"
    r"call\s+dc_lcfo.*?"
    r"solve_dg_overlapping_wannier_generalized_eigenexa.*?"
    r"write_dg_overlapping_wannier_checkpoint",
    main_source,
    re.I | re.S,
), "main_dft needs a concrete construction-to-one-shot-EigenExa-to-checkpoint production adapter"
production_adapter = re.search(
    r"subroutine\s+run_dg_overlapping_wannier_ground_state_for_main(?P<body>.*?)"
    r"end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert production_adapter
adapter_body = production_adapter.group("body")
assert "call build_dg_smooth_partition_of_unity(" in adapter_body.lower(), (
    "production must normalize overlapping core-buffer windows before assembling the global pencil"
)
assert "call assemble_dg_stitched_overlap_density_rows(" in adapter_body.lower(), (
    "production must assemble row-owned overlap and density tiles from normalized buffer coverage"
)
assert "call assemble_dg_stitched_weak_operator_rows(" in source("src/gs/main_dft.f90").lower(), (
    "production Hamiltonian must use partition-gradient-corrected row-owned weak operators"
)
assert "sqrt(ow_partition_weight(p))*ow_box_values(:,p)" in re.sub(
    r"\s+", "", source("src/gs/main_dft.f90").lower()
), "nonlocal projector overlaps must use the same square-root partition weighting"
assert "call symmetrize_dg_distributed_pencil_rows(" in source("src/gs/main_dft.f90").lower(), (
    "production must symmetrize stitched H/S/rho row tiles under the atomic affine action"
)
assert "ow_pencil_translation_subgroup" in source("src/gs/main_dft.f90").lower() and (
    "ow_pencil_coset_representatives" in source("src/gs/main_dft.f90").lower()
), "full affine averaging must use the proven translation-coset factorization"
assert re.search(
    r"prefix\s*=\s*['\"]\./overlapping_wannier_gs['\"]",
    adapter_body,
    re.I,
), (
    "route checkpoint publication must use the communicator-shared run "
    "directory, not a rank-local DC fragment directory"
)
assert not re.search(
    r"prefix\s*=\s*trim\s*\(\s*base_directory\s*\).*?overlapping_wannier_gs",
    adapter_body,
    re.I,
), "route checkpoint prefix must not depend on rank-local base_directory"
assert "checkpoint read rejected:" in adapter_body.lower(), (
    "a rejected route checkpoint must report its exact read/provenance reason"
)
assert not re.search(r"build_dg_core_owned_occupied_subspace\s*\(", adapter_body, re.I), (
    "the production route must not construct a fragment-local occupied direct sum"
)
assert not re.search(
    r"construct_dg_overlapping_wannier_basis\s*\(\s*MPI_COMM_SELF",
    adapter_body,
    re.I,
), "the production route must not construct a fragment-local complement"
assert re.search(
    r"call\s+dc_lcfo\s*\(.*?retained_count\s*=\s*ntarget\s*,\s*&?\s*"
    r"retained_box_count\s*=\s*nstate\s*,\s*&?\s*"
    r"retained_box_contribution\s*=\s*lcfo_fragment_contribution\s*,\s*&?\s*"
    r"retained_occupations\s*=\s*lcfo_retained_occupations",
    adapter_body,
    re.I | re.S,
), "LCFO must retain the finite-temperature window but materialize only occupied buffer rows"
assert re.search(
    r"coefficient_count\s*=\s*max\s*\(\s*dc%nstate_tot\s*,\s*retained_count\s*\)",
    lcfo_source,
    re.I,
), "the in-memory LCFO path must retain requested columns beyond normal nstate_tot"
assert re.search(
    r"subroutine\s+dc_lcfo\s*\(.*?retained_count\s*,\s*"
    r"retained_box_contribution\s*,\s*retained_occupations\s*,\s*write_files\s*,\s*&?\s*"
    r"retained_box_count\s*\)",
    lcfo_source,
    re.I | re.S,
), "the new optional buffer-row count must be appended to preserve positional-call compatibility"
assert re.search(
    r"box_count\s*=\s*retained_count.*?present\s*\(\s*retained_box_count\s*\).*?"
    r"box_count\s*=\s*retained_box_count.*?"
    r"retained_box_contribution\s*\(\s*box_count\s*,\s*product\s*\(\s*nxyz_box\s*\)\s*\)",
    lcfo_source,
    re.I | re.S,
), "LCFO must allocate buffered rows independently of its larger diagonalization window"
assert re.search(
    r"occupations\s*=\s*lcfo_retained_occupations\s*\(\s*1\s*:\s*nstate\s*\)",
    adapter_body,
    re.I,
), "Galerkin occupations must preserve the authoritative LCFO occupation spectrum"
assert re.search(
    r"abs\s*\(\s*sum\s*\(\s*occupations\s*\)\s*-\s*dc%elec_num_tot\s*\)",
    adapter_body,
    re.I,
), "retained LCFO occupations must be rechecked against the total electron count"
assert re.search(
    r"subroutine\s+assign_dg_overlapping_wannier_occupations\b",
    construction_source,
    re.I,
), "retained occupation assignment must live in the overlapping-Wannier namespace"
assert re.search(
    r"maxval\s*\(\s*lcfo_total_symmetry_residual\s*\)\s*>\s*dg_ow_symmetry_tolerance",
    adapter_body,
    re.I,
), "full-affine closure must be accepted before the LCFO seed space is orthonormalized"
composition_call = adapter_body.lower().find(
    "call compose_dg_buffered_orbital_tile_to_physical_grid"
)
fixed_center_call = adapter_body.lower().find("call prepare_ow_fixed_center_group")
affine_measurement_call = adapter_body.lower().find(
    "call measure_dg_rank_fixed_symmetry_residuals_eigenexa"
)
assert 0 <= composition_call < fixed_center_call < affine_measurement_call, (
    "buffer composition must precede fixed-center adaptation and every full-affine proof"
)
translation_adaptation_call = adapter_body.lower().find(
    "call build_dg_group_averaged_occupied_candidates_eigenexa", fixed_center_call
)
point_cogroup_adaptation_call = adapter_body.lower().find(
    "call build_dg_cocycle_averaged_occupied_candidates_eigenexa", translation_adaptation_call + 1
)
post_adaptation_affine_measurement_call = adapter_body.lower().find(
    "call measure_dg_rank_fixed_symmetry_residuals_eigenexa", point_cogroup_adaptation_call
)
assert 0 <= translation_adaptation_call < point_cogroup_adaptation_call < post_adaptation_affine_measurement_call, (
    "occupied adaptation must close translations before the cocycle-aware point cogroup"
)
assert "translation_adapted_occupied" in adapter_body.lower(), (
    "production must publish a separate translation-subgroup adaptation receipt"
)
assert "point_cogroup_adapted_occupied" in adapter_body.lower(), (
    "production must publish a separate cocycle-aware point-cogroup adaptation receipt"
)
assert re.search(
    r"build_dg_cocycle_averaged_occupied_candidates_eigenexa\s*\(.*?"
    r"global_symmetry_map\s*\(\s*:\s*,\s*global_translation_subgroup\s*\).*?"
    r"global_symmetry_map\s*\(\s*:\s*,\s*global_point_representatives\s*\).*?"
    r"global_point_cogroup_product.*?global_translation_cocycle",
    adapter_body,
    re.I | re.S,
), "production point-cogroup adaptation must consume the proven affine factorization"
assert re.search(
    r"nstate\s*>\s*huge\s*\(\s*nstate\s*\)\s*/\s*size\s*\(\s*global_translation_subgroup\s*\)",
    adapter_body,
    re.I,
), "translation-orbit dimension must be overflow-checked before EigenExa initialization"
assert "call accumulate_dg_lcfo_buffer_contributions_to_core" not in adapter_body.lower(), (
    "fragment-core truncation must not define the support of a pre-Wannier symmetry proof"
)
assert "w90_anchors=global_seed_values" in re.sub(r"\s+", "", adapter_body.lower()), (
    "Wannier90 projections must use the complete LCFO core-plus-buffer seed basis"
)
assert re.search(r"call\s+apply_dg_w90_gamma_transform", adapter_body, re.I), (
    "the MLWF transform must be applied in the global LCFO space"
)
assert not re.search(r"allocate\s*\(\s*candidate\s*\(",adapter_body,re.I), (
    "production must not materialize a separate fragment-eigenstate candidate window"
)
assert re.search(
    r"subroutine\s+assemble_dg_distributed_candidate_symmetry",
    construction_source,
    re.I,
), "construction needs streaming distributed candidate symmetry assembly"
assert re.search(
    r"if\s*\(\s*\.not\.\s*distributed_candidates\s*\)\s*then\s*"
    r"call\s+mpi_allgatherv\s*\(\s*candidate_value",
    construction_source,
    re.I | re.S,
), (
    "only the focused replicated-input path may gather candidate boxes"
)
assert "materialize_ow_distributed_core_to_buffer" in adapter_body.lower(), (
    "production must stream the symmetry-closed core into local buffers"
)
assert "materialize_ow_global_tails" not in adapter_body.lower(), (
    "production must not all-gather full-system real-space Wannier tails"
)
assert not re.search(
    r"allocate\s*\(\s*candidate\s*\(\s*ncandidate\s*,\s*nbox",
    adapter_body,
    re.I,
), (
    "production must not zero-pad every fragment candidate onto every rank"
)
assert "prepare_ow_global_point_action" in adapter_body.lower(), (
    "production symmetry must be derived from the full-system crystallographic catalog"
)
assert "prepare_ow_exact_fragment_symmetry" not in adapter_body.lower(), (
    "fragment site symmetry must not be a production prerequisite"
)
canonical_index_body = re.search(
    r"integer\s+function\s+dc_to_canonical_index(?P<body>.*?)end\s+function",
    main_source,
    re.I | re.S,
)
assert canonical_index_body and re.search(
    r"modulo\s*\(\s*index\s*-\s*1\s*,\s*core_count\s*\+\s*2\s*\*\s*buffer_count\s*\)",
    canonical_index_body.group("body"),
    re.I,
), "periodic-buffer projector support must wrap before canonical indexing"
atom_map_body = re.search(
    r"subroutine\s+map_dc_atom_to_physical_atom(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert atom_map_body and re.search(
    r"epsilon\s*\(\s*1d0\s*\).*?maxval\s*\(\s*abs\s*\(\s*dc%system_tot%primitive_a",
    atom_map_body.group("body"),
    re.I | re.S,
), "periodic physical-atom matching needs a scale-relative floating-point tolerance"
assert "pp%nrps_ao" in main_source.lower()
assert "pp%upptbl_ao" in main_source.lower()
assert re.search(
    r"do\s+ll\s*=\s*0\s*,\s*1.*?"
    r"radial_projector.*?pp%upptbl_ao\s*\(\s*1\s*:\s*pp%nrps_ao\s*\(\s*species\s*\)"
    r"\s*,\s*atomic_orbital_ordinals\s*\(\s*ll\s*\+\s*1\s*,\s*species\s*\)\s*,\s*species\s*\)",
    main_source,
    re.I | re.S,
), "every complete s+p seed channel must use the matching PP pseudo-atomic orbital"
projector_adapter_body = re.search(
    r"subroutine\s+build_ow_complete_sp_projectors(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
).group("body")
assert not re.search(r"radial_projector.*?pp%udvtbl", projector_adapter_body, re.I | re.S), (
    "the nonlocal PP projector belongs to the Hamiltonian and must not become an s+p orbital seed"
)
assert "gaussian" not in re.search(
    r"subroutine\s+build_ow_complete_sp_projectors(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
).group("body").lower(), "complete-s+p production seeds must not use a Gaussian fallback"
assert not re.search(r"modulo\s*\(\s*rank\s*\+\s*isym", adapter_body, re.I), (
    "communicator-rank arithmetic is not a physical fragment symmetry"
)
assert "assemble_dg_stitched_weak_operator_rows" in main_source.lower(), (
    "production Hamiltonian must use the boundary-correct stitched weak assembly"
)
hamiltonian_adapter = re.search(
    r"subroutine\s+ow_build_hamiltonian(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert hamiltonian_adapter and not re.search(
    r"call\s+hpsi", hamiltonian_adapter.group("body"), re.I
), "production Hamiltonian must not project the strong fragment stencil"
assert re.search(
    r"if\s*\(\s*yn_dg_dc_overlapping_wannier\s*/=\s*'y'\s*\.and\..*?checkpoint_gs",
    main_source,
    re.I | re.S,
), "overlapping-Wannier route must suppress normal shutdown checkpoint publication"

dispatch_block = re.search(
    rf"if\s*\(\s*{flag}\s*==\s*'y'\s*\)\s*then(?P<body>.*?)else\s+if",
    main_source,
    re.I | re.S,
)
assert dispatch_block
assert not re.search(
    r"\b(?:dc_lcfo|finalize_eigenexa|publish_dg|checkpoint_gs|main_tddft)\b",
    dispatch_block.group("body"),
    re.I,
), "new-route dispatch must not invoke a forbidden stage"

# Disabled behavior is preserved structurally: the conventional DC publication
# ladder remains an else-if ladder and the SCF driver is not globally diverted.
assert re.search(
    r"if\s*\(\s*yn_dc\s*==\s*'y'\s*\)\s*then.*?"
    r"else\s+if\s*\(\s*yn_dc_lcfo_flux\s*==\s*'y'\s*\).*?"
    r"else\s+if\s*\(\s*yn_dc_lcfo\s*==\s*'y'\s*\)",
    main_source,
    re.I | re.S,
), "normal DC LCFO dispatch must remain present"
assert "yn_dg_dc_overlapping_wannier" not in scf_source, (
    "Task 1 must stop before, not alter, the conventional SCF call graph"
)
for workspace in (
    "rho_s_1d",
    "rho_s_sp_1d",
    "exc_1d",
    "eexc_1d",
    "vexc_1d",
    "vexc_sp_1d",
):
    assert re.search(
        rf"if\s*\(\s*allocated\s*\(\s*{workspace}\s*\)\s*\)\s*then"
        rf".*?size\s*\(\s*{workspace}\s*,?\s*1?\s*\)\s*/=\s*nl"
        rf".*?deallocate\s*\(\s*{workspace}\s*\)",
        xc_source,
        re.I | re.S,
    ), f"PZ workspace {workspace} must resize after buffer-local XC evaluation"

assert "dg_overlapping_wannier_types.f90" in dc_cmake
assert "dg_overlapping_wannier_w90.f90" in dc_cmake
assert re.search(r'set\s*\(\s*WANNIER90_COMMS\s+"serial"', w90_builder_source, re.I), (
    "rank-zero Wannier90 library mode requires a serial bundled library"
)
assert re.search(r'set\s*\(\s*WANNIER90_BUILD_TARGETS\s+wannier\s+lib', w90_builder_source, re.I), (
    "bundled serial build must retain both normal-DC executable and OW library"
)
assert not re.search(r'set\s*\(\s*WANNIER90_COMMS\s+"mpi"', w90_builder_source, re.I), (
    "bundled library must not switch to MPI COMMS when SALMON itself uses MPI"
)
for token in (
    "type,public :: s_dg_wannier_tail",
    "type,public :: s_dg_overlapping_wannier_basis",
    "owned_core_physical_ids",
    "physical_grid_ids",
    "gradient",
    "generation",
    "geometry_fingerprint",
    "basis_fingerprint",
    "checked_dg_wannier_extent_product",
    "MPI_Allreduce",
    "MPI_Allgatherv",
):
    assert token.lower() in types_source.lower(), f"missing Task 2 metadata contract: {token}"

assert "dg_overlapping_wannier_construction.f90" in dc_cmake
for forbidden in ("dc_lcfo", "dg_wpw", "direct_sipg"):
    assert forbidden not in construction_source.lower(), (
        f"construction path must not call forbidden stage: {forbidden}"
    )
assert construction_source.lower().count("eigen_pdsyevd_ex_distributed_blocks") == 3, (
    "EigenExa may enter construction only through one import and the two OW-distributed eigensystems"
)
averaged_projector_body = re.search(
    r"subroutine\s+build_dg_group_averaged_occupied_candidates_eigenexa(?P<body>.*?)end\s+subroutine",
    construction_source,
    re.I | re.S,
)
assert averaged_projector_body
averaged_projector_lower = averaged_projector_body.group("body").lower()
assert "orbit_basis" not in averaged_projector_lower, (
    "group averaging must stream occupied images instead of retaining the real-space orbit tensor"
)
assert "cyclic_gram(info%nrow_local,info%ncol_local)" in re.sub(
    r"\s+", "", averaged_projector_lower
), "production orbit Gram must use the EigenExa cyclic distributed layout"
assert re.search(
    r"call\s+zgemm\s*\(\s*'c'\s*,\s*'t'",
    construction_source,
    re.I,
), "periodic symmetry overlaps must use the BLAS matrix product"
assert "candidate_symmetry_raw_unitarity_defect" in construction_source.lower()
assert "candidate_symmetry_polar_correction" in construction_source.lower()
assert "retained_symmetry_projector_gap" in construction_source.lower()
assert "align_dg_fragment_wannier_gauge" in construction_source.lower()
assert "replicate_dg_fragment_wannier_representative" in construction_source.lower()
assert "verify_dg_fragment_wannier_streaming_closure" in construction_source.lower()
assert "build_dg_core_owned_occupied_subspace" in construction_source.lower()
assert "periodic_buffer_boundary_value_norm" in construction_source.lower()
assert "minimum_pivot=" in adapter_body.lower()
assert "stitched overlap lost positive-definite rank" in source(
    "src/gs/dc/dg_overlapping_wannier_metric.f90"
).lower()
assert "ow_collective_operator_fingerprint" in adapter_body.lower()
assert "comm_is_root(nproc_id_global)" not in adapter_body.lower(), (
    "the route communicator root must use its MPI rank, not mutable global rank state"
)
assert re.search(
    r"call\s+write_ow_ground_state_evidence\s*\(\s*spectrum\s*,\s*noccupied\s*,\s*nproc\s*,\s*rank",
    adapter_body,
    re.I,
), "route evidence publication must be owned by rank zero exactly once"
potential_update = re.search(
    r"subroutine\s+dg_dc_update_potential_from_density(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert potential_update
for field in ("dc%rho_tot%f", "dc%Vh_tot%f", "dc%Vxc_tot(is)%f", "dc%vloc_tot(is)%f"):
    assert not re.search(
        rf"all\s*\(\s*ieee_is_finite\s*\(\s*{re.escape(field)}\s*\)\s*\)",
        potential_update.group("body"),
        re.I,
    ), f"distributed DC field {field} must not validate unowned storage"
assert all(token in potential_update.group("body").lower() for token in (
    "dc%mg_tot%is(1)", "dc%mg_tot%ie(1)",
    "dc%mg_tot%is(2)", "dc%mg_tot%ie(2)",
    "dc%mg_tot%is(3)", "dc%mg_tot%ie(3)",
)), "DC potential validation must cover exactly the owned total-system slab"
assert "call build_ow_complete_sp_projectors" in adapter_body.lower()
assert re.search(
    r"global_projection_count\s*=\s*size\s*\(\s*manifest_channels\s*\).*?"
    r"local_target_count\s*=\s*global_projection_count\s*/\s*nproc",
    adapter_body,
    re.I | re.S,
), "global complete-s+p catalog must be rank balanced before LCFO selection"
assert re.search(
    r"ntarget\s*=\s*nstate\s*\+\s*global_projection_count.*?"
    r"call\s+dc_lcfo\s*\(.*?retained_count\s*=\s*ntarget.*?"
    r"retained_box_count\s*=\s*nstate.*?"
    r"retained_box_contribution\s*=\s*lcfo_fragment_contribution",
    adapter_body,
    re.I | re.S,
), "LCFO must separate the finite-temperature occupation window from occupied buffer materialization"
assert re.search(
    r"complete_sp_core_atom_count\s*=\s*dc%system_tot%nion\s*/\s*nproc.*?"
    r"local_target_count\s*/=\s*4\s*\*\s*complete_sp_core_atom_count",
    adapter_body,
    re.I | re.S,
), "production must enforce four complete s+p channels per core-owned atom"
assert "local_seed_overlap" not in adapter_body.lower(), (
    "complete-s+p rank selection must not allocate a replicated dense pre-Wannier localizer"
)
for forbidden in ("augmented_candidate", "augmented_occupied", "occupied_coefficients"):
    assert forbidden not in adapter_body.lower(), (
        f"production must not retain fragment-local complement data: {forbidden}"
    )
assert not re.search(
    r"local_target_count\s*=\s*(?:merge\s*\(\s*)?dg_ow_target_wanniers_per_fragment",
    adapter_body,
    re.I,
), "numeric target input must not override the complete-s+p manifest"
assert not re.search(
    r"noccupied\s*>\s*local_target_count",
    adapter_body,
    re.I,
), "occupied rank must not be pre-cut against the raw shell-channel count"
projector_adapter = re.search(
    r"subroutine\s+build_ow_complete_sp_projectors\b(?P<body>.*?)"
    r"end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert projector_adapter, "missing production complete-s+p pseudopotential adapter"
for token in ("pp%rad", "pp%upptbl_ao", "pp%nrps_ao"):
    assert token.lower() in projector_adapter.group("body").lower(), (
        f"production complete-s+p seeds must use pseudo-atomic orbital data: {token}"
    )
assert "gaussian" not in projection_source.lower(), (
    "periodic complete-s+p projectors must not silently fall back to Gaussian seeds"
)
assert "expected_core_count" in adapter_body.lower() and "ow_partition_weight" in adapter_body.lower(), (
    "stitched metric coverage must be checked against every global physical-grid point"
)
assert re.search(
    r"if\s*\(\s*\.not\.\s*present\s*\(\s*center_representative_box_ids\s*\)\s*\)\s*then.*?"
    r"boundary_value_max\s*>\s*boundary_value_tolerance",
    construction_source,
    re.I | re.S,
), "periodic buffer tails must be measured but not rejected as exterior-zero tails"
assert not re.search(r"call\s+replicate_ow_global_symmetry_orbit", adapter_body, re.I), (
    "production must not copy a representative-fragment gauge across the full system"
)
materialize_position = adapter_body.lower().index("call materialize_ow_distributed_core_to_buffer")
localize_position = adapter_body.lower().index("call run_dg_w90_gamma_library")
metric_position = adapter_body.lower().index("call assemble_dg_stitched_overlap_density_rows")
assert materialize_position < localize_position < metric_position, (
    "localization must use streamed local buffers and finish before metric/SCF publication"
)
assert not re.search(
    r"global_seed_values\s*\(\s*1\s*:\s*nstate.*?=\s*augmented_candidate\s*\(\s*1\s*:\s*nstate",
    adapter_body,
    re.I | re.S,
), "fragment eigenstate indices must not be spliced into fictitious full-system KS seeds"
assert re.search(
    r"compose_dg_buffered_orbital_tile_to_physical_grid\s*\(.*?"
    r"lcfo_fragment_contribution\s*,\s*ow_core_ids\s*,\s*global_seed_values",
    adapter_body,
    re.I | re.S,
), "the occupied block must be composed from LCFO fragment-plus-buffer values"
assert re.search(
    r"do\s+projector_tile_first\s*=.*?"
    r"manifest_channels\s*\(\s*projector_tile_first\s*:\s*projector_tile_last\s*\).*?"
    r"compose_dg_buffered_orbital_tile_to_physical_grid.*?"
    r"global_seed_values\s*\(\s*nstate\s*\+\s*projector_tile_first",
    adapter_body,
    re.I | re.S,
), "the complete-s+p complement must be evaluated and buffer-composed in bounded global-channel tiles"
assert "global_candidate_localizer" not in adapter_body.lower(), (
    "accepted affine-closed seeds must not re-enter the removed dense fixed-rank selector"
)
assert "w90_workspace_peak" in adapter_body and "w90_coordinator_bytes" in adapter_body, (
    "production must publish Wannier90 coordinator and assembly memory receipts"
)
assert "subroutine assemble_dg_periodic_spread_gradient" in localization_source.lower(), (
    "localization must differentiate the same complete periodic spread used by its line search"
)
assert re.search(r"maximum_point_block\s*=\s*4096", localization_source, re.I), (
    "complete spread-gradient assembly must bound its local-grid temporary storage"
)
assert re.search(r"local_action\s*=\s*local_action\s*\+\s*matmul", localization_source, re.I), (
    "complete spread-gradient actions must use blocked dense matrix products"
)
assert not re.search(r"density_potential\s*\(\s*nwannier\s*,\s*npoint\s*\)", localization_source, re.I), (
    "complete spread-gradient assembly must not duplicate the whole local buffer"
)
exponential_body = re.search(
    r"subroutine\s+exponentiate_antihermitian_block\b(?P<body>.*?)end\s+subroutine",
    localization_source,
    re.I | re.S,
)
assert exponential_body and "zheev" in exponential_body.group("body").lower(), (
    "dense global gauge exponential must use one Hermitian eigensolve, not a Taylor matmul series"
)
closure_call = adapter_body.find("call orthonormalize_dg_distributed_seed_space")
localization_call = adapter_body.find("call run_dg_w90_gamma_library")
assert closure_call >= 0, (
    "production OW GS must orthonormalize the accepted full-affine-closed seed space"
)
point_action_body = re.search(
    r"subroutine\s+prepare_ow_global_point_action\b(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert point_action_body and "duplicate_operation" in point_action_body.group("body").lower(), (
    "global point action must deduplicate only identical affine operations"
)
assert "solve_dg_affine_common_fixed_point" in point_action_body.group("body").lower(), (
    "global affine action must report common-center availability without requiring it"
)
assert "build_sawf_operation_index" in point_action_body.group("body").lower(), (
    "full affine product metadata must use the bounded SAWF operation hash index"
)
assert "lookup_sawf_operation_product" in point_action_body.group("body").lower(), (
    "each affine product must be resolved directly from its full normalized key"
)
assert not re.search(
    r"do\s+g\s*=.*?do\s+h\s*=.*?do\s+k\s*=.*?all_maps",
    point_action_body.group("body"),
    re.I | re.S,
), "production must not discover affine products by Nsym-candidate grid-map search"
assert not re.search(r"if\s*\(\s*have_common_center\s*\).*?cycle",point_action_body.group("body"),re.I | re.S), (
    "valid screw/glide operations must not be discarded for lacking a common fixed point"
)
assert localization_call > closure_call, (
    "full-affine closure acceptance and orthonormalization must precede Wannier localization"
)
assert re.search(r"call\s+apply_dg_w90_gamma_transform", adapter_body, re.I), (
    "Wannier90 localization must finish with SALMON canonical Gamma gauge application"
)
assert "replicate_ow_global_symmetry_orbit" not in adapter_body, (
    "production must not replicate a representative fragment after local construction"
)
assert "call run_dg_overlapping_wannier_scf(" not in adapter_body.lower(), (
    "the stitched pencil is a one-shot post-DC solve; repeated OW-SCF would repeatedly "
    "symmetrize H/S/rho and feed its reconstructed density back into the operator"
)
assert adapter_body.lower().count("call solve_dg_overlapping_wannier_generalized_eigenexa(") == 1, (
    "production must solve the symmetrized stitched generalized pencil exactly once with EigenExa"
)
assert "one_shot_operator_fingerprint" in adapter_body.lower(), (
    "one-shot assembly must return a separate fingerprint instead of overwriting the expected operator identity"
)
assert re.search(
    r"one_shot_operator_fingerprint\s*/=\s*operator_fingerprint.*?error\s+stop",
    adapter_body,
    re.I | re.S,
), "one-shot assembly must reject an operator fingerprint mismatch"
adapter_lower = adapter_body.lower()
build_position = adapter_lower.find("call ow_build_hamiltonian(")
solve_position = adapter_lower.find("call solve_dg_overlapping_wannier_generalized_eigenexa(")
density_position = adapter_lower.find("call reconstruct_dg_overlapping_wannier_density(")
checkpoint_position = adapter_lower.find("call write_dg_overlapping_wannier_checkpoint(")
assert 0 <= build_position < solve_position < density_position < checkpoint_position, (
    "production order must be build symmetrized H/S/rho, one-shot EigenExa solve, "
    "density reconstruction, then checkpoint publication"
)
checkpoint_body = re.search(
    r"subroutine\s+populate_ow_checkpoint(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
).group("body")
for stale_refinement_receipt in (
    "refined_residual",
    "refined_orthogonality",
    "refined_condition",
    "refined_coefficients",
    "refined_eigenvalues",
):
    assert stale_refinement_receipt not in checkpoint_body.lower(), (
        f"one-shot checkpoint must not publish stale refinement receipt {stale_refinement_receipt}"
    )
hamiltonian_builder = re.search(
    r"subroutine\s+ow_build_hamiltonian(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert hamiltonian_builder
builder_lower = hamiltonian_builder.group("body").lower()
assemble_position = builder_lower.find("call assemble_dg_stitched_weak_operator_rows(")
symmetrize_position = builder_lower.find("call symmetrize_dg_distributed_pencil_rows(")
assert 0 <= assemble_position < symmetrize_position, (
    "the Hamiltonian builder must symmetrize the assembled stitched pencil before returning it"
)
assert re.search(
    r"inherit_dg_w90_affine_receipts\s*\(.*?w90_closure_defect.*?"
    r"global_retained_group_closure_defect\s*=\s*w90_closure_defect",
    adapter_body,
    re.I | re.S,
), "SCF closure must inherit the accepted full-affine proof through the unitary MLWF gauge"
assert re.search(
    r"call\s+populate_ow_checkpoint\s*\(\s*occupations\s*,\s*condition_number\s*,\s*"
    r"global_retained_group_closure_defect",
    main_source,
    re.I | re.S,
), "V3 checkpoint caller must pass the rank-consistent exact group-algebra closure"
assert re.search(
    r"ow_checkpoint%symmetry_closure_residual\s*=\s*closure_residual",
    main_source,
    re.I,
), "V3 checkpoint writer must publish its explicit closure contract"
for evidence in (
    "localization_initial_spread",
    "localization_final_spread",
    "localization_maximum_gradient",
    "localization_iterations",
    "localization_converged",
):
    assert evidence in ow_checkpoint_source, f"V3 checkpoint misses {evidence}"
assert re.search(
    r"call\s+run_dg_w90_gamma_library.*?if\s*\(\s*\.not\.\s*ok\s*\).*?error\s+stop",
    adapter_body,
    re.I | re.S,
), "rejected Wannier90 MLWF gauges must not reach V3 publication"
assert not re.search(
    r"call\s+align_dg_fragment_wannier_gauge",
    adapter_body,
    re.I,
), "production must not assume independently selected retained spaces differ only by gauge"
assert "inherit_dg_w90_affine_receipts" in adapter_body, (
    "production must preserve the authoritative global action under MLWF gauge rotation"
)
materialize_body = re.search(
    r"subroutine\s+materialize_ow_distributed_core_to_buffer\b(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert materialize_body
assert "mpi_bcast(owner_values" in materialize_body.group("body").lower()
assert "mpi_allgather" not in materialize_body.group("body").lower()
assert "sort_ow_id_positions" in materialize_body.group("body").lower()
assert "find_sorted_ow_id" in materialize_body.group("body").lower()
assert "findloc" not in materialize_body.group("body").lower()
for required in (
    "build_dg_balanced_orbital_ownership",
    "transpose_dg_spatial_cores_to_orbital_owners",
    "redistribute_dg_owned_orbitals_to_center_fragments",
    "assign_dg_periodic_centers_to_fragments",
):
    assert required in construction_source.lower(), f"missing orbital redistribution primitive: {required}"
assert construction_source.lower().count("mpi_alltoallv") >= 2, (
    "both spatial-to-orbital and orbital-to-center-fragment transposes must use MPI_Alltoallv"
)
assert not re.search(r"mod\s*\(\s*ntarget\s*,\s*nproc\s*\)\s*/=\s*0", adapter_body, re.I), (
    "production orbital ownership must support target ranks not divisible by MPI size"
)
assert re.search(
    r"redistribute_dg_owned_orbitals_to_center_fragments\s*\(.*?physical_ids",
    adapter_body,
    re.I | re.S,
), "center-fragment redistribution must send the complete periodic core+buffer box"
assert re.search(
    r"nstate\s*=\s*ceiling\s*\(\s*0\.5d0\s*\*\s*dc%elec_num_tot\s*\)",
    adapter_body,
    re.I,
), "production occupied rank must come from the total electron count in the global LCFO ordering"
assert not re.search(
    r"do\s+source\s*=\s*1\s*,\s*total_box.*?"
    r"do\s+j\s*=\s*1\s*,\s*ncandidate\s*;\s*do\s+i\s*=\s*1\s*,\s*ncandidate.*?"
    r"spatial_overlap\s*\(\s*i\s*,\s*j\s*\)\s*=",
    construction_source,
    re.I | re.S,
), "periodic symmetry overlaps must not use the cubic Fortran scalar loop"

assert "dg_overlapping_wannier_operators.f90" in dc_cmake
for required in ("gradients", "local_potential", "unique-core", "MPI_Allreduce"):
    assert required.lower() in operators_source.lower(), (
        f"missing weak-operator contract: {required}"
    )
for forbidden in ("direct_sipg", "face_hamiltonian", "buffer_volume"):
    assert forbidden not in operators_source.lower(), (
        f"weak operator must not add an independent path: {forbidden}"
    )

for required in (
    "run_dg_overlapping_wannier_scf",
    "unmixed_density_residual",
    "mix_density",
    "rollback_transaction",
):
    assert required.lower() in ow_scf_source.lower(), f"missing Task 8 SCF contract: {required}"
for required in (
    "manifest_magic",
    "shard_magic",
    "versioned_shard_name",
    "call rename",
    "unmixed_density_residual",
    "orthogonality_defect",
    "metric_condition",
    "global_lcfo_fingerprint",
    "occupation_block_fingerprint",
    "affine_cocycle_fingerprint",
    "redistribution_fingerprint",
    "gs_acceptance_receipts",
    "gs_acceptance_tolerance",
):
    assert required.lower() in ow_checkpoint_source.lower(), (
        f"missing Task 8 route-checkpoint contract: {required}"
    )
for forbidden in ("direct_sipg", "lcfo", "eigenexa", "dg_wpw", "checkpoint_gs", "main_tddft"):
    assert forbidden not in ow_scf_source.lower()
    if forbidden != "lcfo":
        assert forbidden not in ow_checkpoint_source.lower()

checkpoint_population = re.search(
    r"subroutine\s+populate_ow_checkpoint\b(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert checkpoint_population
assert not re.search(
    r"MPI_Bcast\s*\(\s*owner_basis\s*,\s*nstate\s*\*\s*nlocal",
    construction_source,
    re.I,
), "production symmetry action must not broadcast every owner's complete basis"
assert not re.search(
    r"allocate\s*\([^)]*lcfo_occupied_representation\s*\(",
    adapter_body,
    re.I | re.S,
), "production must not retain the unused LCFO representation for every point operation"
assert not re.search(
    r"allocate\s*\([^)]*image\s*\(\s*nstate\s*,\s*nlocal\s*\)",
    construction_source,
    re.I | re.S,
), "symmetry residual measurement must tile its image workspace"
assert not re.search(
    r"allocate\s*\([^)]*residual\s*\(\s*nstate\s*,\s*nlocal\s*\)",
    construction_source,
    re.I | re.S,
), "symmetry residual measurement must tile its residual workspace"
direct_eigenexa = re.search(
    r"subroutine\s+eigen_pdsyevd_ex_distributed_blocks\b(?P<body>.*?)end\s+subroutine",
    eigenexa_source,
    re.I | re.S,
)
assert direct_eigenexa, "missing direct distributed-block EigenExa adapter"
assert re.search(r"call\s+eigen_sx\s*\(", direct_eigenexa.group("body"), re.I), (
    "distributed-block EigenExa adapter must diagonalize its local cyclic block directly"
)
assert not re.search(r"\bh\s*\(\s*:\s*,\s*:\s*\)", direct_eigenexa.group("body"), re.I), (
    "distributed-block EigenExa adapter must not accept a replicated dense input"
)
assert "work_matrix" not in direct_eigenexa.group("body").lower(), (
    "distributed-block EigenExa adapter must consume its local block without a duplicate"
)
eigenexa_initializer = re.search(
    r"subroutine\s+init_eigenexa\b(?P<body>.*?)end\s+subroutine",
    source("src/gs/eigenexa_module.f90"),
    re.I | re.S,
)
assert eigenexa_initializer and "direct_block_only" in eigenexa_initializer.group("body").lower(), (
    "OW-sized EigenExa initialization must bypass GS orbital redistribution metadata"
)
assert re.search(
    r"call\s+init_eigenexa_mod\s*\([^\n]*direct_block_only\s*=\s*\.true\.",
    adapter_body,
    re.I,
), "production OW must request direct-block-only EigenExa initialization"
assert re.search(
    r"call\s+init_eigenexa_mod\s*\(\s*info\s*,\s*size\s*\(\s*global_seed_values\s*,\s*1\s*\)\s*,"
    r"[^\n]*direct_block_only\s*=\s*\.true\.",
    adapter_body,
    re.I,
), (
    "the affine residual EigenExa descriptor must use the measured production-seed rank; "
    "initializing it with the smaller occupied rank can return false zero closure/workspace"
)
rank_fixed_residuals = re.search(
    r"subroutine\s+measure_dg_rank_fixed_symmetry_residuals_eigenexa\b(?P<body>.*?)end\s+subroutine",
    construction_source,
    re.I | re.S,
)
assert rank_fixed_residuals, "missing production distributed rank-fixed affine residual implementation"
assert re.search(
    r"assemble_dg_eigenexa_cyclic_metric_block\s*\([^)]*?info%nrow_local\s*,\s*info%ncol_local",
    rank_fixed_residuals.group("body"),
    re.I | re.S,
), (
    "cyclic metric assembly must allocate EigenExa's padded local matrix shape, "
    "not the smaller unpadded cyclic ownership count"
)
for failure_text in (
    "distributed rank-fixed eigensystem",
    "distributed rank-fixed occupied metric is singular",
):
    failure_branch = re.search(
        rf"if\s*\([^\n]*\)\s*then(?P<body>.*?)message\s*=\s*'{re.escape(failure_text)}",
        rank_fixed_residuals.group("body"),
        re.I | re.S,
    )
    assert failure_branch and "ok=.false." in re.sub(r"\s+", "", failure_branch.group("body")).lower(), (
        f"{failure_text} must not inherit a successful metric-assembly status"
    )
    assert "workspace_peak_bytes" in failure_branch.group("body").lower(), (
        f"{failure_text} must publish the workspace measured before rejection"
    )
assert re.search(
    r"call\s+measure_dg_rank_fixed_symmetry_residuals_eigenexa\s*\(",
    adapter_body,
    re.I,
), "production OW construction must use the EigenExa-distributed residual path"
assert "call select_dg_group_generators(" in adapter_body.lower(), (
    "production must prove full affine closure through a deterministic generating set"
)
assert not re.search(
    r"call\s+measure_dg_rank_fixed_symmetry_residuals\s*\(",
    adapter_body,
    re.I,
), "production OW construction must not call the dense-reference residual path"
for replicated_name in (
    "local_metric",
    "metric",
    "metric_vectors",
    "metric_inverse_sqrt",
    "local_overlap",
    "global_overlap",
):
    assert not re.search(
        rf"allocate\s*\([^)]*\b{replicated_name}\s*\(\s*nstate\s*,\s*nstate\s*\)",
        rank_fixed_residuals.group("body"),
        re.I | re.S,
    ), f"Task 3B forbids replicated nstate-square {replicated_name} storage"
assert "assemble_dg_eigenexa_cyclic_metric_block" in rank_fixed_residuals.group("body"), (
    "Task 3B must assemble the rank-fixed metric directly in cyclic EigenExa ownership"
)
assert "eigen_pdsyevd_ex_distributed_blocks" in rank_fixed_residuals.group("body"), (
    "Task 3B must diagonalize the cyclic metric without a replicated dense input"
)
assert not re.search(
    r"MPI_Allreduce\s*\([^\n]*nstate\s*\*\s*nstate",
    rank_fixed_residuals.group("body"),
    re.I,
), "Task 3B forbids all-replicating dense metric or overlap blocks"
symmetry_closed_builder = re.search(
    r"subroutine\s+build_dg_distributed_symmetry_closed_basis\b(?P<body>.*?)end\s+subroutine",
    construction_source,
    re.I | re.S,
)
assert symmetry_closed_builder
assert not re.search(
    r"do\s+left\s*=\s*1\s*,\s*noperation.*?do\s+right\s*=\s*1\s*,\s*noperation.*?"
    r"do\s+source_owner\s*=",
    symmetry_closed_builder.group("body"),
    re.I | re.S,
), "symmetry-closed construction must not repeat an O(Nsym^2*Ngrid) action-table proof"
row_owned_overlap = re.search(
    r"subroutine\s+assemble_dg_distributed_basis_symmetry_overlap_rows\b"
    r"(?P<body>.*?)end\s+subroutine",
    construction_source,
    re.I | re.S,
)
assert row_owned_overlap, "missing row-owned all-operation symmetry overlap assembly"
assert "exchange_dg_point_permuted_orbital_rows" in row_owned_overlap.group("body"), (
    "row-owned symmetry overlaps must use sparse point exchange"
)
assert not re.search(r"MPI_Bcast\s*\(", row_owned_overlap.group("body"), re.I), (
    "row-owned symmetry overlap assembly must not broadcast a complete basis"
)
assert not re.search(
    r"allocate\s*\([^)]*\(\s*nbasis\s*,\s*nbasis\s*,\s*nsym",
    row_owned_overlap.group("body"),
    re.I | re.S,
), "row-owned symmetry overlap assembly must not allocate replicated Nsym*Norb**2"
assert re.search(
    r"subroutine\s+validate_dg_streamed_affine_representation\b",
    construction_source,
    re.I,
), "missing bounded-memory one-operation-at-a-time affine validator"
assert re.search(
    r"subroutine\s+gather_dg_single_symmetry_representation\b",
    construction_source,
    re.I,
), "missing one-operation fixed-center representation gather"
dmn_begin = adapter_body.lower().find("call begin_sawf_dmn(")
dmn_append = adapter_body.lower().find("call append_sawf_dmn_operation(")
dmn_finish = adapter_body.lower().find("call finish_sawf_dmn(")
w90_setup = adapter_body.lower().find("call setup_dg_w90_gamma_library(")
assert min(dmn_begin, dmn_append, dmn_finish, w90_setup) >= 0, (
    "production must publish a streamed fixed-center DMN before Wannier90 setup"
)
assert dmn_begin < dmn_append < dmn_finish < w90_setup, (
    "fixed-center DMN transaction must complete before Wannier90 setup"
)
assert "site_symmetry = .true." in w90_source.lower(), (
    "Wannier90 setup must enable site_symmetry"
)
assert "symmetrize_eps" in w90_source.lower(), (
    "Wannier90 setup must set an explicit strict symmetrize_eps"
)
assert "inquire(file=trim(seed)//'.dmn'" in w90_source.replace(" ", "").lower(), (
    "Wannier90 setup must reject a missing DMN before library entry"
)
assert re.search(r"call\s+validate_dg_w90_convergence_log\s*\(", w90_source, re.I), (
    "Wannier90 run must validate the final convergence receipt"
)
assert re.search(
    r"call\s+gather_dg_single_symmetry_representation\s*\(", adapter_body, re.I
), "production fixed-center DMN must gather one representation operation at a time"
assert "fixed_center_group_order>48" in adapter_body.replace(" ", "").lower(), (
    "production must reject a fixed-center subgroup above crystallographic order 48"
)
assert re.search(r"fixed_center.*inversion", adapter_body, re.I | re.S), (
    "production must require inversion in the fixed-center subgroup"
)
assert re.search(
    r"call\s+inherit_dg_w90_affine_receipts\s*\(",
    adapter_body,
    re.I,
), "production must inherit the accepted full-affine proof through the unitary MLWF gauge"
assert not re.search(
    r"call\s+validate_dg_streamed_affine_representation\s*\(", adapter_body, re.I
), "production must not rebuild all 1536 post-MLWF representation matrices"
assert not re.search(r"w90_symmetry_rows\s*\(", adapter_body, re.I), (
    "production must not retain Nsym row-owned representation matrices"
)
assert not re.search(
    r"assemble_dg_distributed_basis_symmetry_overlap_rows\s*\([^;]*global_symmetry_map",
    adapter_body,
    re.I | re.S,
), "production must not assemble the full-affine all-operation row tensor"
for provenance in (
    "global_lcfo_fingerprint",
    "occupation_block_fingerprint",
    "affine_cocycle_fingerprint",
    "redistribution_fingerprint",
):
    assert re.search(
        rf"ow_checkpoint\s*%\s*{provenance}\s*=",
        checkpoint_population.group("body"),
        re.I,
    ), f"production V3 population omits {provenance}"
assert "ow_checkpoint%gs_acceptance_receipts=" in checkpoint_population.group("body").lower(), (
    "production V3 population omits reconstructed-GS acceptance receipts"
)
assert re.search(
    r"ow_checkpoint\s*%\s*operator_fingerprint\s*=\s*operator_fingerprint\b",
    checkpoint_population.group("body"),
    re.I,
), (
    "route checkpoint provenance must retain the immutable pre-SCF operator "
    "fingerprint, not the final density-dependent Hamiltonian rebuild fingerprint"
)
operator_fingerprint_body = re.search(
    r"function\s+dg_dc_operator_fingerprint\b(?P<body>.*?)end\s+function",
    main_source,
    re.I | re.S,
)
assert operator_fingerprint_body
assert "pp%udvtbl" in operator_fingerprint_body.group("body").lower(), (
    "operator provenance must include the nonlocal pseudopotential operator"
)
assert re.search(
    r"do\s+ii\s*=\s*1\s*,\s*pp\s*%\s*nrps\s*\(\s*kk\s*\).*?"
    r"do\s+jj\s*=\s*0\s*,\s*sum\s*\(\s*pp\s*%\s*nproj\s*\(\s*:\s*,\s*kk\s*\)\s*\)\s*-\s*1",
    operator_fingerprint_body.group("body"),
    re.I | re.S,
), (
    "operator provenance must hash only active nonlocal radial/projector "
    "entries, excluding unused allocated table tails"
)
assert "pp%upptbl_ao" not in operator_fingerprint_body.group("body").lower(), (
    "atomic-orbital projection seeds belong to the basis fingerprint, not "
    "Hamiltonian operator provenance"
)

ow_publication = re.search(
    r"if\s*\(\s*yn_dg_dc_overlapping_wannier\s*==\s*'y'\s*\)\s*then"
    r"(?P<body>.*?)else\s+if",
    main_source,
    re.I | re.S,
)
assert ow_publication
assert re.search(
    r"if\s*\(\s*\.not\.\s*\(\s*sum1\s*<\s*threshold\s*\)\s*\)",
    ow_publication.group("body"),
    re.I,
), (
    "overlapping-Wannier publication must reject an unconverged conventional DC state"
)
assert "run_dg_overlapping_wannier_ground_state_for_main" in ow_publication.group("body")

print("overlapping-Wannier route contract: PASS")
