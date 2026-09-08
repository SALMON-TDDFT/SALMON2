#!/usr/bin/env python3
"""Check the mutable XC payload contract through the density-only wrapper."""

from pathlib import Path
import re


root = Path(__file__).resolve().parents[2]
source = (root / "src/xc/salmon_xc.f90").read_text()
density = re.search(
    r"subroutine exchange_correlation_density\b(.*?)end subroutine exchange_correlation_density",
    source,
    re.I | re.S,
).group(1)
general = re.search(
    r"subroutine exchange_correlation\b(.*?)end subroutine exchange_correlation",
    source,
    re.I | re.S,
).group(1)

intent = r"type\s*\(\s*s_dft_system\s*\)\s*,\s*intent\s*\(\s*inout\s*\)\s*::\s*system"
assert re.search(intent, density, re.I), "density wrapper must expose system payload mutation"
assert re.search(intent, general, re.I), "lower-level XC routine must expose system payload mutation"
assert "system%xc_payload%use_tau_operator = .false." in general

for relative in (
    "src/rt/main_tddft.f90",
    "src/gs/initialization_dft.f90",
    "src/rt/initialization_rt.f90",
    "src/io/write.f90",
    "src/gs/scf_iteration.f90",
    "src/rt/time_evolution_step.f90",
    "src/gs/dc/dcdft.f90",
):
    caller = (root / relative).read_text()
    for match in re.finditer(r"call\s+exchange_correlation(?:_density)?\s*\(([^,\n]+)", caller, re.I):
        actual = match.group(1).strip()
        assert re.fullmatch(r"[a-z_]\w*(?:%[a-z_]\w*)*", actual, re.I), (
            f"{relative}: mutable system actual is not a definable variable: {actual}"
        )

print("PASS XC mutable payload intent and definable caller contract")
