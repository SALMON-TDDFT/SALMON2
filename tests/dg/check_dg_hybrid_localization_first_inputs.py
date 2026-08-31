#!/usr/bin/env python3
"""Contract for localization-first Hybrid energy and DC-seed inputs."""

import argparse
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
GLOBAL = (ROOT / "src/io/salmon_global.f90").read_text(errors="replace").lower()
INPUT = (ROOT / "src/io/inputoutput.f90").read_text(errors="replace").lower()
SI64 = (
    ROOT
    / "tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_continuation.in"
).read_text(errors="replace").lower()


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


for declaration in (
    r"character\(1\)\s*::\s*yn_dg_hybrid_continuation_scf",
    r"real\(8\)\s*::\s*dg_hybrid_symmetry_energy_window",
    r"character\(16\)\s*::\s*dg_dc_seed_mode",
    r"character\(256\)\s*::\s*dg_dc_seed_directory",
):
    assert re.search(declaration, GLOBAL), f"missing global declaration: {declaration}"

dc_namelist = INPUT[INPUT.index("namelist/dc/") : INPUT.index("!! == default for &unit")]
for name in (
    "yn_dg_hybrid_continuation_scf",
    "dg_hybrid_symmetry_energy_window",
    "dg_dc_seed_mode",
    "dg_dc_seed_directory",
):
    assert name in dc_namelist, f"{name} must belong to &dc"

for default in (
    "yn_dg_hybrid_continuation_scf = 'n'",
    "dg_hybrid_symmetry_energy_window = -1d0",
    "dg_dc_seed_mode = 'off'",
    "dg_dc_seed_directory = ''",
):
    assert default in INPUT, f"missing compatibility default: {default}"

input_broadcasts = INPUT[
    INPUT.index("!! == bcast for dc") : INPUT.index("!! == bcast for dg_fragment")
]
for name in (
    "yn_dg_hybrid_continuation_scf",
    "dg_hybrid_symmetry_energy_window",
    "dg_dc_seed_mode",
    "dg_dc_seed_directory",
):
    assert f"call comm_bcast({name}" in input_broadcasts, (
        f"missing global broadcast for {name}"
    )

dc_log_start = INPUT.index("'dc', inml_dc")
dc_log = INPUT[dc_log_start : INPUT.index("close(fh_variables_log)", dc_log_start)]
for name in (
    "yn_dg_hybrid_continuation_scf",
    "dg_hybrid_symmetry_energy_window",
    "dg_dc_seed_mode",
    "dg_dc_seed_directory",
):
    assert f"'{name}'" in dc_log or f'"{name}"' in dc_log, (
        f"missing variables.log entry for {name}"
    )

packed = compact(INPUT).replace("&", "")
window_conversion = (
    "if(dg_hybrid_symmetry_energy_window>=0d0)"
    "dg_hybrid_symmetry_energy_window="
    "dg_hybrid_symmetry_energy_window*uenergy_to_au"
)
assert window_conversion in packed, (
    "the nonnegative energy window must be converted independently while preserving -1d0"
)
assert "ieee_is_finite(dg_hybrid_symmetry_energy_window)" in packed
assert "dg_hybrid_symmetry_energy_window/=-1d0" in packed
assert "dg_hybrid_symmetry_energy_window<0d0" in packed

assert "selectcase(trim(dg_dc_seed_mode))" in packed
assert "case('off','write','read','auto')" in packed
assert "len_trim(dg_dc_seed_directory)==0" in packed
assert "call yn_argument_check(yn_dg_hybrid_continuation_scf)" in INPUT

continuation_check = INPUT[INPUT.index("subroutine check_bad_input"):]
cutoff_branch_start = "if(yn_dg_hybrid_continuation_scf=='y')then"
assert cutoff_branch_start in compact(continuation_check)
continuation_check = compact(continuation_check)
continuation_check = continuation_check[continuation_check.index(cutoff_branch_start) :]
continuation_check = continuation_check.split("endif", 1)[0]
for token in (
    "ieee_is_finite(energy_cut)",
    "ieee_is_finite(lambda_cut)",
    "lambda_cut<=0d0",
    "ieee_is_finite(wannier_pw_cutoff)",
    "wannier_pw_cutoff<=0d0",
):
    assert token in continuation_check, (
        f"Hybrid continuation omits independent cutoff validation: {token}"
    )

for forbidden in (
    "dg_hybrid_symmetry_energy_window=energy_cut",
    "dg_hybrid_symmetry_energy_window=wannier_pw_cutoff",
    "wannier_pw_cutoff=dg_hybrid_symmetry_energy_window",
):
    assert forbidden not in packed, f"independent cutoff controls were aliased: {forbidden}"

match = re.search(r"dg_hybrid_symmetry_energy_window\s*=\s*([^\s!/,]+)", SI64)
assert match, "Si64 production continuation input must set an explicit energy window"
assert float(match.group(1).replace("d", "e")) == 0.2, (
    "Si64 production continuation must start with the reviewed 0.2 a.u. window"
)


def check_invalid_input_exit(executable: Path) -> None:
    cases = (
        (
            "dg_hybrid_symmetry_energy_window=-2d0",
            "dg_hybrid_symmetry_energy_window must be -1 or nonnegative",
        ),
        (
            "dg_hybrid_symmetry_energy_window=NaN",
            "dg_hybrid_symmetry_energy_window must be -1 or nonnegative",
        ),
        (
            "dg_dc_seed_mode='invalid'",
            "dg_dc_seed_mode must be off, write, read, or auto",
        ),
        (
            "dg_dc_seed_mode='read', dg_dc_seed_directory=' '",
            "dg_dc_seed_directory is required when dg_dc_seed_mode is enabled",
        ),
        (
            "yn_dg_hybrid_continuation_scf='y', energy_cut=NaN, "
            "lambda_cut=1d0, wannier_pw_cutoff=1d0",
            "DG continuation requires finite energy_cut",
        ),
        (
            "yn_dg_hybrid_continuation_scf='y', energy_cut=0d0, "
            "lambda_cut=0d0, wannier_pw_cutoff=1d0",
            "DG continuation requires positive finite lambda_cut",
        ),
        (
            "yn_dg_hybrid_continuation_scf='y', energy_cut=0d0, "
            "lambda_cut=NaN, wannier_pw_cutoff=1d0",
            "DG continuation requires positive finite lambda_cut",
        ),
        (
            "yn_dg_hybrid_continuation_scf='y', energy_cut=0d0, "
            "lambda_cut=1d0, wannier_pw_cutoff=NaN",
            "DG continuation requires positive finite wannier_pw_cutoff",
        ),
        (
            "yn_dg_hybrid_continuation_scf='y', energy_cut=0d0, "
            "lambda_cut=1d0, wannier_pw_cutoff=0d0",
            "DG continuation requires positive finite wannier_pw_cutoff",
        ),
        (
            "yn_dg_hybrid_continuation_scf='y', energy_cut=-1d0, "
            "lambda_cut=1d0, wannier_pw_cutoff=1d0",
            "DG continuation requires yn_dg_dc_overlapping_wannier='y'",
        ),
    )
    for setting, expected_message in cases:
        with tempfile.TemporaryDirectory(prefix="salmon-hybrid-input-") as temp:
            result = subprocess.run(
                [str(executable)],
                cwd=temp,
                input=f"&dc\n {setting}\n/\n",
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=30,
            )
            assert result.returncode != 0, (
                f"invalid Hybrid input exited zero: {setting}\n{result.stdout}"
            )
            assert expected_message in result.stdout, (
                f"invalid Hybrid input missed diagnostic: {setting}\n{result.stdout}"
            )


def check_window_unit_conversion(executable: Path) -> None:
    cases = (
        ("27.211386245988d0", 1.0),
        ("-1d0", -1.0),
    )
    for input_value, expected_au in cases:
        with tempfile.TemporaryDirectory(prefix="salmon-hybrid-units-") as temp:
            result = subprocess.run(
                [str(executable)],
                cwd=temp,
                input=(
                    "&units\n unit_system='A_eV_fs'\n/\n"
                    "&dc\n"
                    f" dg_hybrid_symmetry_energy_window={input_value},\n"
                    " dg_dc_seed_mode='invalid'\n/\n"
                ),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=30,
            )
            assert result.returncode != 0
            assert "dg_dc_seed_mode must be off, write, read, or auto" in result.stdout
            variables_log = (Path(temp) / "variables.log").read_text()
            match = re.search(
                r"dg_hybrid_symmetry_energy_window\s*=\s*([-+0-9.e]+)",
                variables_log,
                re.IGNORECASE,
            )
            assert match, "variables.log omits the converted Hybrid energy window"
            assert abs(float(match.group(1)) - expected_au) < 1e-6, (
                "Hybrid energy-window conversion changed its value or the -1 sentinel"
            )


parser = argparse.ArgumentParser()
parser.add_argument("--salmon-executable", type=Path)
args = parser.parse_args()
if args.salmon_executable:
    executable = args.salmon_executable.resolve()
    check_invalid_input_exit(executable)
    check_window_unit_conversion(executable)

print("PASS localization-first Hybrid input contracts")
