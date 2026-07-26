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
    r"call\s+dispatch_dg_overlapping_wannier_route.*?"
    r"(?:return|error\s+stop)",
    main_source,
    re.I | re.S,
), "main_dft needs an explicit, terminating new-route dispatch"
assert re.search(
    r"subroutine\s+dispatch_dg_overlapping_wannier_route.*?"
    r"error\s+stop\s*['\"]overlapping Wannier route: construction not implemented['\"]",
    dcdft_source,
    re.I | re.S,
), "Task 1 dispatcher must stop before forbidden stages"

dispatch_block = re.search(
    rf"if\s*\(\s*{flag}\s*==\s*'y'\s*\)\s*then(?P<body>.*?)end\s*if",
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

print("overlapping-Wannier route contract: PASS")
