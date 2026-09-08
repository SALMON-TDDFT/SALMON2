#!/usr/bin/env python3
"""Production contract for SCDM-seeded fragment-WF checkpoint reuse."""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
GLOBAL = (ROOT / "src/io/salmon_global.f90").read_text().lower()
INPUT = (ROOT / "src/io/inputoutput.f90").read_text().lower()
MAIN = (ROOT / "src/gs/main_dft.f90").read_text().lower()
BUILDER = (ROOT / "src/gs/dc/dg_hybrid_fragment_wannier.f90").read_text().lower()
ADMISSION = (ROOT / "src/gs/dc/dg_hybrid_fragment_admission.f90").read_text().lower()
FIXTURE = (ROOT / "tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in").read_text().lower()
SMOKE_RUNNER_PATH = ROOT / "tests/dg/run_dg_fragment_wf_production_smoke.py"

for declaration in (
    r"character\(16\)\s*::\s*dg_fragment_wf_checkpoint_mode",
    r"character\(256\)\s*::\s*dg_fragment_wf_checkpoint_directory",
):
    assert re.search(declaration, GLOBAL), f"missing control declaration: {declaration}"

dc_namelist = INPUT[INPUT.index("namelist/dc/"):INPUT.index("!! == default for &unit")]
for name in ("dg_fragment_wf_checkpoint_mode", "dg_fragment_wf_checkpoint_directory"):
    assert name in dc_namelist
    assert f"call comm_bcast({name}" in INPUT
    assert f"'{name}'" in INPUT or f'"{name}"' in INPUT
assert "dg_fragment_wf_checkpoint_mode = 'auto'" in INPUT
assert "dg_fragment_wf_checkpoint_directory = 'dg-fragment-wf-checkpoint'" in INPUT
packed_input = re.sub(r"\s+|&", "", INPUT)
assert "selectcase(trim(dg_fragment_wf_checkpoint_mode))" in packed_input
assert "case('off','write','read','auto')" in packed_input
assert "case('scdm','spectral','random')" in packed_input
assert "dg_fragment_w90_initial_projection='scdm'" in packed_input

assert "use dg_fragment_scdm_gauge" in BUILDER
assert "use dg_fragment_wf_checkpoint" in BUILDER
construct = BUILDER.split("subroutine construct_fragment_wannier", 1)[1].split(
    "end subroutine construct_fragment_wannier", 1
)[0]
assert construct.index("build_dg_fragment_scdm_gauge") < construct.index("setup_dg_w90_gamma_library")
assert "precomputed_a_matrix=initial_a_matrix" in re.sub(r"\s+|&", "", construct)

build_match = re.search(
    r"subroutine\s+build_dg_hybrid_fragment_wannier\s*\(.*?"
    r"end\s+subroutine\s+build_dg_hybrid_fragment_wannier\b",
    BUILDER,
    re.S,
)
assert build_match
build = build_match.group(0)
for call in (
    "probe_dg_fragment_wf_checkpoint",
    "decide_dg_fragment_wf_restart",
    "read_dg_fragment_wf_checkpoint",
    "write_dg_fragment_wf_checkpoint",
):
    assert f"call {call}" in build, f"missing collective restart action: {call}"
assert build.index("probe_dg_fragment_wf_checkpoint") < build.index("construct_fragment_wannier")
assert build.index("read_dg_fragment_wf_checkpoint") < build.index("construct_fragment_wannier")
assert build.index("construct_fragment_wannier") < build.index("write_dg_fragment_wf_checkpoint")
assert "restore_fragment_cache_from_checkpoint" in build
assert "pack_fragment_cache_for_checkpoint" in build
assert "if(trim(wf_policy)/='off')build_required=.true." in re.sub(r"\s+|&", "", build)
assert "wf_policy_min==wf_policy_max" in re.sub(r"\s+|&", "", build)
assert "projection_min==projection_max" in re.sub(r"\s+|&", "", build)
assert '"-attempt-"' in build

route = MAIN.split("subroutine run_dg_hybrid_divided_ground_state_for_main", 1)[1].split(
    "end subroutine run_dg_hybrid_divided_ground_state_for_main", 1
)[0]
call = re.search(r"call\s+build_dg_hybrid_fragment_wannier_from_dc_seed\s*\((.*?)\)", route, re.S)
assert call
for token in (
    "dg_fragment_wf_checkpoint_mode",
    "dg_fragment_wf_checkpoint_directory",
    "dg_fragment_w90_initial_projection",
    "dg_dc_seed_publication_id",
    "dg_dc_seed_contract%ownership_fingerprint",
    "dg_dc_seed_contract%immutable_fingerprint",
):
    assert token in call.group(1), f"production call omits restart identity: {token}"
assert "atomic_create_directory(trim(dg_fragment_wf_checkpoint_directory)" in route
assert "[dg-fragment-wf]" in route
assert "checkpoint_hit=" in route
assert "[dg-hybrid-divided] projected_basis_fingerprint=" in route
terminal_log = route.split("fixed-density/non-self-consistent divided wf+pw lcfo solved once", 1)[1]
assert "electron_defect=" in terminal_log[:800]

restore_start = BUILDER.index("subroutine restore_fragment_cache_from_checkpoint")
restore_end = BUILDER.index("end subroutine restore_fragment_cache_from_checkpoint", restore_start)
restore = BUILDER[restore_start:restore_end]
assert "contract%rank/=rank" not in restore, (
    "a total-communicator checkpoint rank must not be compared with the rank in a one-rank fragment communicator"
)

# Keep every DC seed visible until the energy-aware initializer has sorted the
# spectrum and completed a cutoff or degenerate boundary.  Pre-slicing by
# initial_count+guard_count can silently split a symmetry multiplet.
compact_admission = re.sub(r"\s+|&", "", ADMISSION)
assert "trial_seed_count" not in compact_admission
assert "coefficients,reference%energies,initial_count,guard_count" in compact_admission

for setting in (
    "dg_dc_seed_mode='auto'",
    "dg_dc_seed_directory='../dc-seed'",
    "dg_fragment_wf_checkpoint_mode='auto'",
    "dg_fragment_wf_checkpoint_directory='../fragment-wf-checkpoint'",
    "dg_fragment_w90_initial_projection='scdm'",
):
    assert setting in re.sub(r"\s+", "", FIXTURE)

assert SMOKE_RUNNER_PATH.is_file(), "missing short Si8 production restart smoke runner"
smoke_runner = SMOKE_RUNNER_PATH.read_text().lower()
for token in (
    "si8_overlapping_wannier",
    'not hit["checkpoint_hit"]',
    'miss["checkpoint_hit"]',
    "projected_basis_fingerprint",
    "electron_defect",
    "dg_fragment_wf.manifest",
    "incomplete",
    "--seed-directory",
    "occupied_checkpoint_fingerprint",
    'incomplete["wannier_fragment_ids"] != expected_fragment_ids',
    "struct.unpack_from",
    "seed_preloaded",
    "--analyze-existing",
    "--evidence-output",
    "require_mpi_completion",
    "parse_finite_float",
    "compare_smoke_runs",
):
    assert token in smoke_runner, f"short production smoke omits {token}"

print("PASS production SCDM and fragment-WF restart route")
