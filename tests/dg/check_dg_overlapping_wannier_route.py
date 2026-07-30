#!/usr/bin/env python3
"""Source contract for the isolated overlapping-Wannier DC route."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def source(path: str) -> str:
    return (ROOT / path).read_text()


global_source = source("src/io/salmon_global.f90")
input_source = source("src/io/inputoutput.f90")
main_source = source("src/gs/main_dft.f90")
scf_source = source("src/gs/scf_iteration_dft.f90")
dcdft_source = source("src/gs/dc/dcdft.f90")
types_source = source("src/gs/dc/dg_overlapping_wannier_types.f90")
construction_source = source("src/gs/dc/dg_overlapping_wannier_construction.f90")
projection_source = source("src/gs/dc/dg_overlapping_wannier_projection.f90")
operators_source = source("src/gs/dc/dg_overlapping_wannier_operators.f90")
ow_scf_source = source("src/gs/dc/dg_overlapping_wannier_scf.f90")
ow_solver_source = source("src/gs/dc/dg_overlapping_wannier_solver.f90")
ow_checkpoint_source = source("src/gs/dc/dg_overlapping_wannier_checkpoint.f90")
dc_cmake = source("src/gs/dc/CMakeLists.txt")
xc_source = source("src/xc/salmon_xc.f90")
si64_runner_source = source("tests/dg/run_si64_overlapping_wannier_gate.py")
si64_checker_source = source("tests/dg/check_si64_overlapping_wannier_gate.py")

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
):
    assert tolerance in global_source
    assert tolerance in input_source
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
    (r"\byn_eigenexa\s*==\s*'y'", "EigenExa forbidden"),
    (r"\byn_dg_wpw_production\s*==\s*'y'", "WPW forbidden"),
    (r"\byn_dg_dc_local_periodic\s*==\s*'y'", "direct-DG forbidden"),
    (r"\byn_dg_fragment_rt\s*==\s*'y'", "conventional RT forbidden"),
    (r"\byn_self_checkpoint\s*==\s*'y'", "normal checkpoint forbidden"),
    (r"\bcheckpoint_interval\s*>=\s*1", "periodic checkpoint forbidden"),
]
for condition, requirement in route_checks:
    assert re.search(
        rf"if\s*\(\s*{flag}\s*==\s*'y'.*?{condition}",
        input_source,
        re.I | re.S,
    ), f"overlapping-Wannier validation must enforce: {requirement}"

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
    r"construct_dg_overlapping_wannier_basis.*?"
    r"run_dg_overlapping_wannier_scf.*?"
    r"write_dg_overlapping_wannier_checkpoint",
    main_source,
    re.I | re.S,
), "main_dft needs a concrete construction-to-SCF-to-checkpoint production adapter"
production_adapter = re.search(
    r"subroutine\s+run_dg_overlapping_wannier_ground_state_for_main(?P<body>.*?)"
    r"end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert production_adapter
adapter_body = production_adapter.group("body")
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
assert re.search(
    r"build_dg_core_owned_occupied_subspace\s*\(.*?system%rocc",
    adapter_body,
    re.I | re.S,
), "Wannier occupied inclusion must use the DC core-weighted chemical-potential subspace"
assert re.search(
    r"assign_dg_dc_local_basis_occupations\s*\(\s*dc%elec_num_tot\s*,\s*occupations",
    adapter_body,
    re.I,
), "Galerkin occupations must use the DC total core electron count"
assert re.search(
    r"allocate\s*\(\s*candidate\s*\(\s*local_candidate_count\s*,\s*nbox",
    adapter_body,
    re.I,
), (
    "production candidates must remain fragment-local before tail materialization"
)
assert re.search(
    r"call\s+construct_dg_overlapping_wannier_basis\s*\(\s*mpi_comm_self",
    adapter_body,
    re.I | re.S,
), "production must close Wannier construction inside each periodic buffer fragment"
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
assert "materialize_ow_global_tails" in adapter_body.lower(), (
    "production must gather only retained fragment Wannier tails"
)
assert not re.search(
    r"allocate\s*\(\s*candidate\s*\(\s*ncandidate\s*,\s*nbox",
    adapter_body,
    re.I,
), (
    "production must not zero-pad every fragment candidate onto every rank"
)
assert "build_dc_translation_symmetry_map" in adapter_body.lower(), (
    "production symmetry must be derived from DC fragment geometry"
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
assert "assemble_dg_overlapping_wannier_weak_operators" in main_source.lower(), (
    "production Hamiltonian must use the Task 5 weak unique-core assembly"
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
    r"if\s*\(\s*yn_dg_dc_local_periodic\s*/=\s*'y'\s*\.and\.\s*"
    r"yn_dg_dc_overlapping_wannier\s*/=\s*'y'\s*\.and\..*?checkpoint_gs",
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
for forbidden in ("dc_lcfo", "eigenexa", "dg_wpw", "direct_sipg"):
    assert forbidden not in construction_source.lower(), (
        f"construction path must not call forbidden stage: {forbidden}"
    )
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
assert "global_metric_rejected_rank" in adapter_body.lower()
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
    r"local_target_count\s*=\s*size\s*\(\s*manifest_channels\s*\)",
    adapter_body,
    re.I,
), "complete-s+p manifest count must be recorded before direct-sum selection"
assert re.search(
    r"construct_dg_overlapping_wannier_basis\s*\(.*?"
    r"projection_seed_values\s*=\s*manifest_values.*?"
    r"local_target_count\s*=\s*ow_basis\s*%\s*target_rank",
    adapter_body,
    re.I | re.S,
), "measured occupied plus complete-shell direct-sum rank must become the production target"
assert re.search(
    r"complete_sp_core_atom_count\s*=\s*count\s*\(\s*manifest_channels\s*%\s*l\s*==\s*0\s*\).*?"
    r"local_target_count\s*/=\s*4\s*\*\s*complete_sp_core_atom_count",
    adapter_body,
    re.I | re.S,
), "production must enforce four complete s+p channels per core-owned atom"
assert re.search(
    r"construct_dg_overlapping_wannier_basis\s*\(.*?"
    r"projection_seed_values\s*=\s*manifest_values",
    adapter_body,
    re.I | re.S,
), "production target selection must consume the complete-s+p projector values"
assert re.search(
    r"augmented_candidate\s*\(\s*1\s*:\s*local_candidate_count\s*,\s*:\s*\)"
    r"\s*=\s*candidate.*?"
    r"augmented_candidate\s*\(\s*local_candidate_count\s*\+\s*1\s*:"
    r"\s*local_candidate_count\s*\+\s*local_target_count\s*,\s*:\s*\)"
    r"\s*=\s*manifest_values",
    adapter_body,
    re.I | re.S,
), "raw complete-s+p orbitals must be appended to the fragment eigenstate candidates"
assert re.search(
    r"periodic_box_gradients\s*\(\s*augmented_candidate.*?"
    r"augmented_gradient",
    adapter_body,
    re.I | re.S,
), "the appended raw s+p orbitals need buffer-periodic finite-difference gradients"
assert re.search(
    r"augmented_occupied\s*=\s*\(\s*0d0\s*,\s*0d0\s*\).*?"
    r"augmented_occupied\s*\(\s*1\s*:\s*local_candidate_count\s*,\s*:\s*\)"
    r"\s*=\s*occupied_coefficients",
    adapter_body,
    re.I | re.S,
), "occupied coefficients must be zero-extended over the appended s+p shell"
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
assert re.search(
    r"allocate\s*\(\s*pairs\s*\(\s*ntarget\s*,\s*ncore\s*\)\s*\)",
    adapter_body,
    re.I,
), "metric ownership mask must cover every local unique-core quadrature point"
assert re.search(
    r"if\s*\(\s*\.not\.\s*present\s*\(\s*center_representative_box_ids\s*\)\s*\)\s*then.*?"
    r"boundary_value_max\s*>\s*boundary_value_tolerance",
    construction_source,
    re.I | re.S,
), "periodic buffer tails must be measured but not rejected as exterior-zero tails"
assert re.search(
    r"call\s+replicate_dg_fragment_wannier_representative",
    adapter_body,
    re.I,
), "production must materialize every translated fragment from one symmetry representative"
assert not re.search(
    r"call\s+align_dg_fragment_wannier_gauge",
    adapter_body,
    re.I,
), "production must not assume independently selected retained spaces differ only by gauge"
assert re.search(
    r"call\s+verify_dg_fragment_wannier_streaming_closure",
    adapter_body,
    re.I,
), "production must verify authoritative retained tails with bounded streaming"
assert re.search(
    r"MPI_Allgather\s*\(\s*dc\s*%\s*i_frag.*?rank_fragment.*?"
    r"center_owner_fragment.*?rank_fragment\s*\(\s*source\s*\)",
    main_source,
    re.I | re.S,
), "materialized Wannier fragment IDs must be independent of MPI rank order"
materialize_body = re.search(
    r"subroutine\s+materialize_ow_global_tails\b(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert materialize_body
assert "source_position" in materialize_body.group("body").lower()
assert "sort_ow_id_positions" in materialize_body.group("body").lower()
assert "findloc" not in materialize_body.group("body").lower()
assert not re.search(
    r"do\s+point\s*=\s*1\s*,\s*nbox.*?count\s*\(",
    materialize_body.group("body"),
    re.I | re.S,
), "tail materialization must not perform quadratic physical-ID lookup"
assert re.search(
    r"call\s+build_dg_core_owned_occupied_subspace",
    adapter_body,
    re.I,
), "production occupied rank must come from the DC core-weighted subspace"
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
):
    assert required.lower() in ow_checkpoint_source.lower(), (
        f"missing Task 8 route-checkpoint contract: {required}"
    )
for forbidden in ("direct_sipg", "lcfo", "eigenexa", "dg_wpw", "checkpoint_gs", "main_tddft"):
    assert forbidden not in ow_scf_source.lower()
    assert forbidden not in ow_checkpoint_source.lower()

checkpoint_population = re.search(
    r"subroutine\s+populate_ow_checkpoint\b(?P<body>.*?)end\s+subroutine",
    main_source,
    re.I | re.S,
)
assert checkpoint_population
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
