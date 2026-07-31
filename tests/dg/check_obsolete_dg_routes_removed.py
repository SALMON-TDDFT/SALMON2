#!/usr/bin/env python3
"""Source contract for removing obsolete experimental DG execution routes."""

from __future__ import annotations

import fnmatch
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SELF = Path(__file__).resolve()

FORBIDDEN_INPUTS = {
    "yn_adaptive_basis",
    "yn_adaptive_basis_dg",
    "yn_adaptive_basis_update",
    "yn_dg_dc_local_periodic",
    "yn_dg_expdiag_delta_h",
    "yn_dg_expdiag_global_field",
    "yn_dg_expdiag_global_flux",
    "yn_dg_expdiag_project_h",
    "yn_dg_expdiag_refresh_fixed_func",
    "yn_dg_expdiag_xi_split",
    "yn_dg_frag",
    "yn_dg_fragment_from_dcdft",
    "yn_dg_fragment_rt",
    "yn_dg_full_h_eigen_seed",
    "yn_dg_full_h_wannier_band_gauge",
    "yn_dg_hse_ace",
    "yn_dg_lcfo_seed_exhaustive_check",
    "yn_dg_mixed_z",
    "yn_dg_mixed_z_decomp_output",
    "yn_dg_mixed_z_include_ww",
    "yn_dg_mixed_z_local_current_writeback_total",
    "yn_dg_mixed_z_local_prop_writeback",
    "yn_dg_mixed_z_local_pz_writeback_total",
    "yn_dg_mixed_z_local_rho_writeback_wwonly",
    "yn_dg_nodal_rt",
    "yn_dg_subspace_diag",
    "yn_dg_wpw_checkpoint_rt",
    "yn_dg_wpw_fixed_h_relaxation",
    "yn_dg_wpw_global_projected_correction",
    "yn_dg_wpw_h_epsilon_s_correction",
    "yn_dg_wpw_metric_preconditioner",
    "yn_dg_wpw_preconditioner",
    "yn_dg_wpw_production",
    "yn_dg_wpw_s_orthogonal_pw",
    "yn_dg_wpw_search_history",
    "basis_update_threshold",
    "dg_bpw_auto",
    "dg_bpw_auto_accuracy",
    "dg_bpw_auto_max_n",
    "dg_bpw_auto_min_n",
    "dg_bpw_auto_report",
    "dg_bpw_position_mode",
    "dg_mixed_z_local_prop_backend",
    "dg_mixed_z_neighbor_env_shell",
    "dg_mixed_z_frag_local_field_block",
    "dg_mixed_z_direct_origin",
    "dg_mixed_z_ww_position_branch",
    "dg_mixed_z_wp_position_branch",
    "dg_mixed_z_pp_position_branch",
    "dg_mixed_z_polarization_branch",
    "dg_nodal_gs_relax_step",
    "dg_nodal_gs_max_iter",
    "dg_nodal_gs_tol",
    "dg_nodal_taylor_order",
    "dg_subspace_extra_states",
    "dg_subspace_pw_vectors",
    "dg_wpw_checkpoint_manifest",
    "dg_wpw_checkpoint_rank_prefix",
    "dg_wpw_checkpoint_identity_tolerance",
    "dg_wpw_exp_max_corrector",
    "dg_wpw_exp_corrector_tolerance",
    "dg_wpw_exp_norm_tolerance",
    "dg_wpw_extra_states",
    "dg_wpw_gap_threshold",
    "dg_wpw_global_correction_max_iterations",
    "dg_wpw_global_correction_restart",
    "dg_wpw_global_correction_state_batch",
    "dg_wpw_global_correction_tolerance",
    "dg_wpw_metric_cutoff",
    "dg_wpw_scf_max_iter",
    "dg_wpw_scf_mix",
    "dg_wpw_scf_residual_tolerance",
    "dg_wpw_window_buffer",
    "dg_wpw_window_width",
    "eps_dg_frag",
    "k_cutoff_plane_wave",
    "n_plane_waves_dg",
    "niter_dg_frag_rt_max",
    "time_integrator_dg_fragment",
    "yn_plane_wave_basis",
}

FORBIDDEN_SOURCE_PREFIXES = (
    "src/common/dg_wpw_",
    "src/gs/dc/dg_wpw_",
    "src/rt/dg/rt_dg_fragment",
    "src/rt/dg/rt_dg_integrator_",
    "src/rt/dg/rt_dg_wpw_",
    "src/rt/dg/rt_dg_nodal_",
)

FORBIDDEN_SOURCE_FILES = {
    "src/common/dg_generalized_algebra.f90",
    "src/common/dg_dc_direct_sipg.f90",
    "src/common/dg_ground_state_checkpoint.f90",
    "src/poisson/poisson_dg_distributed.f90",
    "src/gs/dc/dg_dc_ground_state.f90",
    "src/gs/dc/dg_dc_ground_state_adapter.f90",
    "src/gs/dc/dg_dc_handoff.f90",
    "src/gs/dc/dg_dc_local_basis_ground_state.f90",
    "src/gs/dc/dg_wannier_pw_scf.f90",
}

FORBIDDEN_ASSET_GLOBS = (
    "tests/dg/check_direct_amn_global_gauge.py",
    "tests/dg/check_full_h_*",
    "tests/dg/check_global_wannier_generalized_eigen.py",
    "tests/dg/check_krylov_velocity_integrator.py",
    "tests/dg/check_mixed_z_field_sign.py",
    "tests/dg/check_momentum_active_axis_openmp.py",
    "tests/dg/check_velocity_gauge_full_h_seed.py",
    "tests/dg/*dg_dc_direct*",
    "tests/dg/*dg_dc_ground_state*",
    "tests/dg/*dg_dc_handoff*",
    "tests/dg/*dg_dc_local_basis*",
    "tests/dg/*dg_dc_local_periodic*",
    "tests/dg/*dg_ground_state_checkpoint*",
    "tests/dg/*dg_wpw*",
    "tests/dg/*nodal*",
    "tests/dg/*wpw*",
    "tests/dg/fixtures/*wpw*",
    "samples/exercise_dg_fragment_C2H2/**",
    "samples/exercise_dg_fragment_rt/**",
    "samples/exercise_dg_rt_hse_test/**",
)

REQUIRED_SOURCES = {
    "src/gs/dc/dg_overlapping_wannier_checkpoint.f90",
    "src/gs/dc/dg_overlapping_wannier_construction.f90",
    "src/rt/dg/rt_dg_overlapping_wannier.f90",
}

REQUIRED_INPUTS = {
    "yn_dg_dc_overlapping_wannier",
    "yn_dg_length_gauge",
    "yn_dg_overlapping_wannier_rt",
    "yn_dg_overlapping_wannier_rt_restart",
}

PROTECTED_NORMAL_SOURCES = {
    "src/gs/dc/lcfo.f90": "src/gs/dc/CMakeLists.txt",
    "src/gs/eigen_subdiag_eigenexa.f90": "src/gs/CMakeLists.txt",
}

SCAN_ROOTS = (
    ROOT / "src",
    ROOT / "tests" / "dg",
    ROOT / "samples",
)

HISTORICAL_NOTICE = "historical/removed"
HISTORICAL_DOC_TOKENS = (
    "samples/exercise_dg_fragment_c2h2",
    "samples/exercise_dg_fragment_rt",
    "samples/exercise_dg_rt_hse_test",
    "yn_dg_frag",
    "yn_dg_fragment_rt",
    "yn_dg_mixed_z",
    "yn_dg_wpw_",
    "rt_dg_fragment",
    "dg_wpw_",
)
CURRENT_REMOVAL_DOCS = {
    "docs/plans/2026-07-31-obsolete-dg-route-inventory.md",
    "docs/plans/2026-07-31-remove-obsolete-experimental-dg-routes-design.md",
    "docs/plans/2026-07-31-remove-obsolete-experimental-dg-routes.md",
    "docs/plans/2026-07-31-remove-obsolete-experimental-dg-routes-results.md",
}
FORBIDDEN_TEST_REGISTRATION_TOKENS = (
    "check_direct_amn_global_gauge",
    "check_full_h_",
    "check_global_wannier_generalized_eigen",
    "check_krylov_velocity_integrator",
    "check_mixed_z_field_sign",
    "check_momentum_active_axis_openmp",
    "check_velocity_gauge_full_h_seed",
)


def repository_files() -> list[Path]:
    files: list[Path] = []
    for scan_root in SCAN_ROOTS:
        files.extend(path for path in scan_root.rglob("*") if path.is_file())
    files.extend(ROOT.glob("**/CMakeLists.txt"))
    return sorted(set(files))


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def text_files() -> list[tuple[Path, str]]:
    result: list[tuple[Path, str]] = []
    for path in repository_files():
        if path.resolve() == SELF:
            continue
        try:
            result.append((path, path.read_text(errors="strict")))
        except (UnicodeDecodeError, OSError):
            continue
    return result


def main() -> int:
    failures: list[str] = []
    all_files = repository_files()
    texts = text_files()
    global_source = (ROOT / "src/io/salmon_global.f90").read_text().lower()
    input_source = (ROOT / "src/io/inputoutput.f90").read_text().lower()
    gs_dispatch_source = (ROOT / "src/gs/main_dft.f90").read_text().lower()
    rt_dispatch_source = (ROOT / "src/rt/main_tddft.f90").read_text().lower()

    for required_input in sorted(REQUIRED_INPUTS):
        if required_input not in global_source:
            failures.append(f"missing retained input declaration: {required_input}")
        if required_input not in input_source:
            failures.append(f"missing retained input handling: {required_input}")

    if "namelist/dg_fragment/" in input_source.replace(" ", ""):
        failures.append("forbidden namelist remains: dg_fragment")
    if "time_evolution_dg_fragment" in rt_dispatch_source:
        failures.append("forbidden RT dispatch remains: time_evolution_dg_fragment")
    if rt_dispatch_source.find(
        "call run_dg_overlapping_wannier_coefficient_rt"
    ) > rt_dispatch_source.find(
        "call initialization_rt"
    ):
        failures.append("coefficient RT dispatch must precede conventional RT initialization")
    if "yn_dg_dc_local_periodic" in gs_dispatch_source or "yn_dg_wpw_production" in gs_dispatch_source:
        failures.append("forbidden GS dispatch remains")

    for required in sorted(REQUIRED_SOURCES):
        if not (ROOT / required).is_file():
            failures.append(f"missing retained source: {required}")

    for source_path, cmake_path in sorted(PROTECTED_NORMAL_SOURCES.items()):
        if not (ROOT / source_path).is_file():
            failures.append(f"missing protected normal source: {source_path}")
            continue
        cmake = (ROOT / cmake_path).read_text()
        if Path(source_path).name not in cmake:
            failures.append(
                f"missing protected normal CMake entry: {cmake_path}:"
                f"{Path(source_path).name}"
            )

    for path in all_files:
        rel = relative(path)
        if rel in FORBIDDEN_SOURCE_FILES:
            failures.append(f"forbidden source file: {rel}")
        if any(rel.startswith(prefix) for prefix in FORBIDDEN_SOURCE_PREFIXES):
            failures.append(f"forbidden source prefix: {rel}")
        if any(fnmatch.fnmatch(rel, pattern) for pattern in FORBIDDEN_ASSET_GLOBS):
            failures.append(f"forbidden focused asset: {rel}")

    for path, text in texts:
        lower = text.lower()
        rel = relative(path)
        if path.name == "CMakeLists.txt":
            for token in FORBIDDEN_TEST_REGISTRATION_TOKENS:
                if token in lower:
                    failures.append(f"forbidden test registration: {rel}: {token}")
        for token in sorted(FORBIDDEN_INPUTS):
            token_pattern = re.compile(rf"(?<![a-z0-9_]){re.escape(token)}(?![a-z0-9_])")
            for match in token_pattern.finditer(lower):
                line = lower.count("\n", 0, match.start()) + 1
                failures.append(f"forbidden input: {rel}:{line}: {token}")

    for docs_root in (ROOT / "docs/plans", ROOT / "docs/superpowers"):
        for path in docs_root.rglob("*.md"):
            rel = relative(path)
            if rel in CURRENT_REMOVAL_DOCS:
                continue
            lower = path.read_text(errors="replace").lower()
            if any(token in lower for token in HISTORICAL_DOC_TOKENS):
                if HISTORICAL_NOTICE not in lower:
                    failures.append(
                        "obsolete executable documentation lacks historical notice: "
                        f"{rel}"
                    )

    if failures:
        print("obsolete DG route removal contract: FAIL")
        for finding in sorted(set(failures)):
            print(f"  {finding}")
        print(f"total forbidden inventory entries: {len(set(failures))}")
        return 1

    print("obsolete DG route removal contract: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
