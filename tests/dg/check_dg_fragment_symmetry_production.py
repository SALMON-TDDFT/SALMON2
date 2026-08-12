#!/usr/bin/env python3
"""Production wiring contract for exact buffered-fragment crystallographic symmetry."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
MAIN = (ROOT / "src/gs/main_dft.f90").read_text().lower()
SYMMETRY = (ROOT / "src/gs/dc/dg_overlapping_wannier_symmetry.f90").read_text().lower()
CONSTRUCTION = (ROOT / "src/gs/dc/dg_overlapping_wannier_construction.f90").read_text().lower()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


for routine in (
    "load_sawf_crystallographic_catalog_auto",
    "build_dg_fragment_site_stabilizer",
    "build_dg_fragment_group_representation",
    "promote_dg_exact_global_subgroup",
    "project_dg_fragment_covariant_operators",
):
    require(re.search(rf"\bcall\s+{routine}\b", MAIN) is not None,
            f"production overlapping-Wannier route does not call {routine}")

require("fingerprint_dg_exact_fragment_symmetry" in MAIN,
        "V3 basis provenance must bind each exact fragment group")
require("exact_fragment_symmetry_fingerprints" in MAIN and "mpi_allgather" in MAIN,
        "V3 basis provenance must bind the rank-ordered fragment-group collection")

require("dc%system_tot%rion" in MAIN and "dc%system_tot%kion" in MAIN,
        "fragment symmetry discovery must use instantaneous coordinates and species")
require(re.search(r"fragment.*atom.*(mask|index)", MAIN) is not None,
        "production must construct a separate buffered-fragment atom set")
construction_call = re.search(
    r"call\s+construct_dg_overlapping_wannier_basis\s*\((.*?)\)\s*\n",
    MAIN, re.S)
require((construction_call is None and "call run_dg_w90_gamma_library" in MAIN and
         "call localize_dg_occupation_blocks" not in MAIN) or
        (construction_call is not None and "local_symmetry_map" not in construction_call.group(1)),
        "fragment-local Wannier generation must not impose a fragment-only point group")
require("point_group_symbol" in MAIN and "space_group_number" in MAIN,
        "production diagnostics must identify each exact fragment group")
require("pre_projection_defect" in MAIN and "post_projection_defect" in MAIN,
        "production evidence must report covariance before and after projection")
require("cross_block_scalar_residual" in MAIN and "cross_block_vector_residual" in MAIN,
        "global promotion must audit scalar and vector cross-fragment blocks")
require("global_inversion_promoted" in MAIN,
        "V3 publication must report whether exact global inversion was promoted")
require("project_ow_exact_global_group" in MAIN,
        "fragment-spanning symmetries must be projected as one exact global group")
require("assemble_dg_distributed_basis_symmetry_overlap" in MAIN,
        "global checkpoint representation must be measured from the actual distributed Wannier gauge")
require("assemble_dg_distributed_basis_symmetry_overlap_rows" in MAIN and
        "validate_dg_factored_point_cogroup_gauge" in MAIN,
        "post-MLWF symmetry validation must remain row-owned and memory bounded")
exact_global_body = re.search(
    r"subroutine\s+project_ow_exact_global_group(.*?)end\s+subroutine", MAIN, re.S)
require(exact_global_body is not None and
        "matmul(metric_inverse,symmetry_overlap" in exact_global_body.group(1),
        "global Wannier action must solve S D = <g phi|phi> before group synchronization")
require("mpi_allgather(ow_core_ids" in exact_global_body.group(1) and
        "call build_dg_pointwise_affine_owner_map" in exact_global_body.group(1) and
        "findloc(all_physical_ids" in CONSTRUCTION,
        "full-system operations must resolve every mapped global grid point to its owner")
require("validate_sawf_fragment_symmetry_map" not in exact_global_body.group(1),
        "fragment permutation must not be required by an exact full-system operation")
translation_builder = re.search(
    r"subroutine\s+build_ow_fragment_permutation_representation(.*?)end\s+subroutine",
    MAIN, re.S)
require(translation_builder is not None and
        translation_builder.group(1).index("call mpi_allgather(target_fragment_local") <
        translation_builder.group(1).index("allocate(product_table"),
        "translation product table must be built after its global fragment map is collected")
require("global_exact_group_promoted" in MAIN and "global_exact_group_order" in MAIN,
        "V3 publication must report the simultaneously promoted exact global group")
require("selected_operations" in MAIN and "fixed_residual" in MAIN,
        "global projection must use every exact operation fixing one inversion center")
require("replicate_ow_global_symmetry_orbit" in MAIN,
        "representative local Wanniers must be propagated by full-system affine symmetry")
require("build_dg_fragment_symmetry_orbits" in SYMMETRY,
        "fragment-local generation must support multiple full-system symmetry orbits")
require("call build_dg_fragment_symmetry_orbits" in MAIN,
        "production replication must select one representative per fragment orbit")
require(MAIN.index("call promote_dg_exact_global_subgroup") <
        MAIN.index("call project_dg_fragment_covariant_operators"),
        "local fragment symmetry must not project global matrices before promotion")
require(MAIN.index("call project_dg_fragment_covariant_operators") <
        MAIN.index("call compute_dg_overlapping_wannier_matrix_fingerprints"),
        "S/H/X/V must be projected before V3 matrix fingerprints are published")

for forbidden in ("thermal_symmetry_tolerance", "parent_point_group", "force_parent_symmetry"):
    require(forbidden not in MAIN and forbidden not in SYMMETRY,
            f"production must not restore thermally broken parent symmetry: {forbidden}")

print("PASS exact buffered-fragment symmetry is wired into production publication")
