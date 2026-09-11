#!/usr/bin/env python3
"""Task 8 production gate; helpers alone do not establish main-route wiring.

This gate intentionally remains RED until the divided production route is
separated from legacy complete-system construction and wired to the new kernels.
It is a source contract, not a replacement for MPI/physical regression tests.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/main_dft.f90").read_text().lower()
source = re.sub(r"!.*", "", source).replace("&", "")
entry = "run_dg_hybrid_divided_ground_state_for_main"
assert re.search(r"if\s*\(\s*yn_dg_hybrid_divided_scf\s*==\s*'y'\s*\)\s*then\s*call\s+" + entry, source), (
    "Task 8 incomplete: divided main must dispatch to its own DC-seed construction entry"
)
assert "bare overlapping-wannier gs is retired" in source
match = re.search(r"\bsubroutine\s+" + entry + r"\b(.*?)\bend subroutine\s+" + entry, source, re.S)
assert match, "missing separated divided production routine"
route = match.group(1)
fragment_builder = (ROOT / "src/gs/dc/dg_hybrid_fragment_wannier.f90").read_text().lower()
admission_source = (ROOT / "src/gs/dc/dg_hybrid_fragment_admission.f90").read_text().lower()
assert "'dgfw'" in route and "base_directory)//'dg_hybrid_fragment_w90'" not in route, (
    "production W90 seed path must remain below the library fixed-length seed limit"
)
assert "'/w'" in fragment_builder and "'/construction_wannier'" not in fragment_builder, (
    "fragment W90 leaf name exceeds the library fixed-length seed budget"
)
certify = fragment_builder.split("subroutine certify_seed_span", 1)[1].split(
    "end subroutine certify_seed_span", 1
)[0]
compress = fragment_builder.split("subroutine metric_compress_candidates", 1)[1].split(
    "end subroutine metric_compress_candidates", 1
)[0]
for expression in (
    "coefficients=matmul",
    "reconstructed=matmul",
):
    assert expression in re.sub(r"\s+", "", certify), (
        "full DC span certification must use the matrix backend, not scalar cubic loops"
    )
for expression in ("gram=matmul", "retained=matmul", "retained_gram=matmul"):
    assert expression in re.sub(r"\s+", "", compress), (
        "fragment metric compression must use the matrix backend"
    )
compact_admission = re.sub(r"\s+|&", "", admission_source)
assert "trial_seed_count" not in compact_admission, (
    "pre-slicing DC seeds can split a degenerate energy boundary"
)
assert "coefficients,reference%energies,initial_count,guard_count" in compact_admission, (
    "the complete seed inventory must reach the energy-aware initializer"
)
assert "nproc==dc%n_frag" in re.sub(r"\s+", "", route), (
    "divided production entry must require MPI size equal to fragment count"
)
compact_route = re.sub(r"\s+", "", route)
assert "fragment_lattice(axis,axis)=dc%system_tot%hgs(axis)*real(raw_grid(axis),8)" in compact_route
assert "fragment_reciprocal_lattice(axis,axis)=2d0*acos(-1d0)/fragment_lattice(axis,axis)" in compact_route
assert not re.search(r"build_dg_hybrid_fragment_wannier_from_dc_seed\(.*?system%primitive_a",
                     route, re.S), "full-system lattice was passed to fragment Wannier90"
required = (
    "build_dg_hybrid_fragment_wannier_from_dc_seed",
    "select_dg_hybrid_core_wannier",
    "prepare_dg_hybrid_selected_catalog",
    "build_dg_hybrid_projected_local_fragment_basis",
    "prepare_dg_hybrid_selected_trial",
    "export_dg_hybrid_selected_basis_frame",
    "prepare_dg_hybrid_schwarz_candidate_inventory",
    "initialize_dg_hybrid_schwarz_state",
    "build_dg_hybrid_schwarz_schedule",
    "initialize_dg_hybrid_interface_continuation",
    "solve_dg_hybrid_generalized_once_and_publish",
)
for name in required:
    assert re.search(r"\bcall\s+" + name + r"\b", route), f"missing production call: {name}"
assert route.index("call build_dg_hybrid_fragment_wannier_from_dc_seed") < route.index(
    "call initialize_dg_hybrid_interface_continuation"
)
ordered = (
    "call build_dg_hybrid_fragment_wannier_from_dc_seed",
    "call select_dg_hybrid_core_wannier",
    "call prepare_dg_hybrid_selected_catalog",
    "call build_dg_hybrid_projected_local_fragment_basis",
    "call prepare_dg_hybrid_selected_trial",
    "call initialize_dg_hybrid_interface_continuation",
)
for before, after in zip(ordered, ordered[1:]):
    assert route.index(before) < route.index(after), f"production order must keep {before} before {after}"
diagnostic_begin = route.index("if(density_diagnostic)then")
diagnostic_end = route.index("call initialize_dg_hybrid_interface_continuation", diagnostic_begin)
diagnostic = route[diagnostic_begin:diagnostic_end]
assert diagnostic.count("call dc_lcfo(") == 1
assert "retained_core_density=diagnostic_core" in diagnostic
assert "write_files=.false." in diagnostic
assert "call dc_lcfo(" not in route[:diagnostic_begin] + route[diagnostic_end:], (
    "ordinary LCFO is allowed only in the optional fixed-density comparison"
)
for name in ("run_dg_overlapping_wannier_ground_state_for_main",
             "setup_dg_w90_gamma_library", "run_dg_w90_gamma_library",
             "apply_dg_hybrid_divided_fragment_hpsi", "solve_dg_hybrid_fragment_spectrum"):
    assert not re.search(r"\bcall\s+" + name + r"\b", route), f"legacy production fallback: {name}"
assert not re.search(r"\b384\b", route), "material-specific state count in production"
for forbidden in (
    r"nstate\s*=\s*[^\n]*retained_rank",
    r"state_count\s*=\s*[^\n]*retained_rank",
    r"coefficients[^\n]*=\s*0[^\n]*pw",
    r"dc_seed_coefficients[^\n]*\(\s*:\s*selected_count",
):
    assert not re.search(forbidden, route), f"raw/unprojected state shortcut remains: {forbidden}"
solver_name = "solve_dg_hybrid_schwarz_fragments"
solver_match = re.search(r"\bsubroutine\s+" + solver_name + r"\b(.*?)\bend subroutine\s+" + solver_name,
                         source, re.S)
assert solver_match, "missing Schwarz fragment solve callback"
solver = solver_match.group(1)
assert solver_name in route, "divided SCF does not receive the Schwarz fragment callback"
assert "advance_dg_hybrid_schwarz_epoch" in solver
assert "assign_dg_hybrid_schwarz_occupations" in solver
assert "dg_hybrid_fragment_cg_steps" in solver, "Schwarz production route ignores the user CG cap"
assert "300d0" in solver, "common-mu occupation epoch omits the 300 K electron temperature"
for obsolete in (
    "solve_dg_hybrid_bounded_fragments",
    "solve_dg_hybrid_divided_fragments",
    "apply_dg_hybrid_divided_fragment_hpsi",
    "apply_dg_hybrid_divided_fragment_metric",
    "assemble_dg_hybrid_divided_core_density",
):
    assert not re.search(r"\bsubroutine\s+" + obsolete + r"\b", source), (
        f"obsolete divided production callback remains: {obsolete}"
    )
print("fragment-local DC-to-Wannier production route: PASS")
