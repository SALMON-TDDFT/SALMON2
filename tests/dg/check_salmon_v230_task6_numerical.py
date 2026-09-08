#!/usr/bin/env python3
"""Focused static contracts for the SALMON v2.3.0/DG numerical merge."""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[2]
TASK6_PATHS = (
    "src/atom/pp/prep_pp.f90",
    "src/atom/pp/salmon_pp.f90",
    "src/common/density_matrix.f90",
    "src/common/hamiltonian.f90",
    "src/common/initialization.f90",
    "src/common/structures.f90",
    "src/common/total_energy.f90",
    "src/gs/gram_schmidt_orth.f90",
    "src/gs/scf_iteration.f90",
    "src/gs/dc/lcfo.f90",
    "src/io/checkpoint_restart.f90",
    "src/io/write.f90",
    "src/poisson/hartree.f90",
    "src/poisson/poisson_periodic.f90",
    "src/rt/initialization_rt.f90",
    "src/rt/time_evolution_step.f90",
    "src/xc/builtin_tbmbj.f90",
    "src/xc/salmon_xc.f90",
)


def source(relative):
    return (ROOT / relative).read_text(encoding="utf-8", errors="replace")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def require_tokens(relative, *tokens):
    text = source(relative).casefold()
    missing = [token for token in tokens if token.casefold() not in text]
    require(not missing, f"{relative}: missing {', '.join(missing)}")


def require_no_markers(relative):
    text = source(relative)
    require(
        not re.search(r"^(?:<<<<<<<|\|\|\|\|\|\|\||=======|>>>>>>>)", text, re.M),
        f"{relative}: unresolved conflict marker",
    )


def hamiltonian():
    relative = "src/common/hamiltonian.f90"
    require_no_markers(relative)
    require_tokens(relative, "apply_nonlocal", "nvtxStartRange('nonlocal potential'")
    text = source(relative)
    require(
        re.search(r"if\s*\(\s*yn_jm\s*==\s*'n'\s*\.and\.\s*apply_nonlocal\s*\)\s*then", text, re.I),
        f"{relative}: nonlocal application is not controlled by apply_nonlocal",
    )


def auto_merged_kernels():
    """Protect the nine clean auto-merges that still combine shared behavior."""
    expectations = {
        "src/atom/pp/prep_pp.f90": ("use nvtx_wrapper", "ppg%uv_so /= ppg%uv_so"),
        "src/atom/pp/salmon_pp.f90": ("use nvtx_wrapper", "rho_nlcc_tmp", "if(flag_cuboid)"),
        "src/common/density_matrix.f90": (
            "use nvtx_wrapper",
            "calc_current_plusU",
            "micro_current",
            "system%rmatrix_B",
        ),
        "src/common/total_energy.f90": ("use nvtx_wrapper", "xc_payload%use_tau_operator", "yn_fix_func"),
        "src/io/write.f90": ("MPI_File_write_all", "write_dg_polarization_data", "write_dg_polarization_response_3d"),
        "src/poisson/poisson_periodic.f90": (
            "use nvtx_wrapper",
            "set_poisson_contract_context",
            "poisson_ft_hse_sr",
            "poisson_ffte_hse_sr",
        ),
        "src/rt/initialization_rt.f90": (
            "use nvtx_wrapper",
            "init_dm_unfold",
            "initialization_rt_common",
            "yn_conventional_from_dcdft",
        ),
        "src/xc/builtin_tbmbj.f90": ("rho_thr=1d-10", "bracketed"),
    }
    for relative, tokens in expectations.items():
        require_no_markers(relative)
        require_tokens(relative, *tokens)


def initialization():
    for relative in ("src/common/initialization.f90", "src/common/structures.f90"):
        require_no_markers(relative)
    require_tokens(
        "src/common/initialization.f90",
        "if_stop = .false.",
        "non-orthogonal lattice",
        "flag_blacs_gridinit",
    )
    require_tokens("src/common/structures.f90", "flag_blacs_gridinit", "if_divide_rspace")


def gs_scf():
    for relative in ("src/gs/gram_schmidt_orth.f90", "src/gs/scf_iteration.f90"):
        require_no_markers(relative)
    require_tokens("src/gs/gram_schmidt_orth.f90", "cublasZdotc", "ZDOTC", "ZDSCAL")
    require_tokens(
        "src/gs/scf_iteration.f90",
        "yn_dc",
        "yn_hse",
        "calc_xc_hse_fft",
        "ace_update_decision",
    )
    text = source("src/gs/scf_iteration.f90")
    require(text.count("call update_vlocal(") == 1, "src/gs/scf_iteration.f90: duplicate update_vlocal call")


def lcfo():
    relative = "src/gs/dc/lcfo.f90"
    require_no_markers(relative)
    require_tokens(
        relative,
        "lcfo_eigensolver",
        "retained_count",
        "retained_box_contribution",
        "retained_occupations",
        "retained_eigenvalues",
        "build_retained_occupations",
        "build_retained_box_contribution",
        "diag_chefsi_driver",
    )
    text = source(relative)
    require(text.count("select case(trim(lcfo_eigensolver))") == 1, f"{relative}: solver selection must occur once")
    for name, guard in (("diag_eigenexa", "USE_EIGENEXA"), ("diag_chefsi_driver", "USE_SCALAPACK")):
        pattern = rf"#ifdef\s+{guard}.*?call\s+{name}\b"
        require(re.search(pattern, text, re.I | re.S), f"{relative}: {name} lacks {guard} guard")


def checkpoint():
    relative = "src/io/checkpoint_restart.f90"
    require_no_markers(relative)
    require_tokens(relative, "use nvtx_wrapper", "yn_reset_occupation_restart", "write_rho0", "read_rho0")


def hartree():
    for relative in ("src/poisson/hartree.f90", "src/poisson/poisson_periodic.f90"):
        require_no_markers(relative)
    require_tokens("src/poisson/hartree.f90", "use nvtx_wrapper", "hse_omega", "SALMON_HSE_SR_HARTREE")
    require_tokens("src/poisson/poisson_periodic.f90", "poisson_periodic", "rho")


def rt():
    for relative in ("src/rt/initialization_rt.f90", "src/rt/time_evolution_step.f90"):
        require_no_markers(relative)
    initialization_rt = source("src/rt/initialization_rt.f90")
    require(
        re.search(r"subroutine\s+initialization_rt\s*\([^)]*ppn\s*,\s*unfold\s*\)", initialization_rt, re.I | re.S),
        "src/rt/initialization_rt.f90: conventional wrapper lost the v2.3 unfold output",
    )
    require(
        re.search(r"call\s+initialization_rt_common\s*\([^)]*ppn\s*,\s*unfold\s*\)", initialization_rt, re.I | re.S),
        "src/rt/initialization_rt.f90: conventional wrapper does not forward unfold",
    )
    require(
        re.search(r"subroutine\s+initialization_rt_dg_hybrid\b.*?type\s*\(\s*s_unfold\s*\)\s*::\s*unfold.*?unfold\s*=\s*unfold", initialization_rt, re.I | re.S),
        "src/rt/initialization_rt.f90: DG wrapper does not supply local unfold state",
    )
    require_tokens(
        "src/rt/time_evolution_step.f90",
        "use nvtx_wrapper",
        "dm_unfold",
        "calc_density_matrix_and_energy_plusU",
        "ace_update_decision",
        "compute_stage2_residual",
    )


def xc():
    for relative in ("src/xc/builtin_tbmbj.f90", "src/xc/salmon_xc.f90"):
        require_no_markers(relative)
    require_tokens(
        "src/xc/salmon_xc.f90",
        "salmon_xctype_r2scan",
        "build_gauge_invariant_tau",
        "xc_func%use_kinetic_energy .or. xc_func%use_current",
        "orbital-dependent exchange-correlation requires orbitals",
        "use nvtx_wrapper",
        "ieee_is_finite",
    )
    require_tokens("src/xc/builtin_tbmbj.f90", "builtin_tbmbj", "tau")


CHECKS = {
    "auto-merged": auto_merged_kernels,
    "hamiltonian": hamiltonian,
    "initialization": initialization,
    "gs-scf": gs_scf,
    "lcfo": lcfo,
    "checkpoint": checkpoint,
    "hartree": hartree,
    "rt": rt,
    "xc": xc,
}


requested = sys.argv[1:]
check_all_paths = not requested
if requested:
    unknown = [name for name in requested if name not in CHECKS]
    require(not unknown, f"unknown check(s): {', '.join(unknown)}")
else:
    requested = list(CHECKS)

for name in requested:
    CHECKS[name]()
    print(f"Task 6 {name}: PASS")

if check_all_paths:
    for relative in TASK6_PATHS:
        require_no_markers(relative)
print("Task 6 numerical integration contracts: PASS")
