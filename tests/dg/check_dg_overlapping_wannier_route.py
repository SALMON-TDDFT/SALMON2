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
operators_source = source("src/gs/dc/dg_overlapping_wannier_operators.f90")
ow_scf_source = source("src/gs/dc/dg_overlapping_wannier_scf.f90")
ow_checkpoint_source = source("src/gs/dc/dg_overlapping_wannier_checkpoint.f90")
dc_cmake = source("src/gs/dc/CMakeLists.txt")

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
    r"ncandidate\s*=\s*(?:system%no|local_candidate_count)\s*\*\s*nproc",
    adapter_body,
    re.I,
), (
    "production candidates must retain the per-fragment direct-sum identity"
)
assert re.search(r"candidate_offset\s*=\s*.*(?:fragment|i_frag)", adapter_body, re.I), (
    "production candidates need an explicit fragment/local-orbital block offset"
)
assert "build_dc_translation_symmetry_map" in adapter_body.lower(), (
    "production symmetry must be derived from DC fragment geometry"
)
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

print("overlapping-Wannier route contract: PASS")
