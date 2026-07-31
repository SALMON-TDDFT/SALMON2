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
    "src/gs/dc/dg_dc_ground_state.f90",
    "src/gs/dc/dg_dc_ground_state_adapter.f90",
    "src/gs/dc/dg_dc_handoff.f90",
    "src/gs/dc/dg_dc_local_basis_ground_state.f90",
    "src/gs/dc/dg_wannier_pw_scf.f90",
}

FORBIDDEN_ASSET_GLOBS = (
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
        for token in sorted(FORBIDDEN_INPUTS):
            token_pattern = re.compile(rf"(?<![a-z0-9_]){re.escape(token)}(?![a-z0-9_])")
            for match in token_pattern.finditer(lower):
                line = lower.count("\n", 0, match.start()) + 1
                failures.append(f"forbidden input: {rel}:{line}: {token}")

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
