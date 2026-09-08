#!/usr/bin/env python3
"""Require normalized exact validation of both W90 projection controls."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/io/inputoutput.f90").read_text(errors="replace").lower()


def violations(source: str) -> list[str]:
    compact = re.sub(r"\s+|&", "", source)
    problems: list[str] = []
    lower = compact.find("callstring_lowercase(dg_ow_w90_initial_projection)")
    validate = compact.find("selectcase(trim(dg_ow_w90_initial_projection))")
    cases = compact.find("case('spectral','random')", validate)
    fatal = compact.find("dg_ow_w90_initial_projectionmustbespectralorrandom", validate)
    if min(lower, validate, cases, fatal) < 0 or not lower < validate < cases < fatal:
        problems.append("global W90 projection must be normalized then validated as spectral|random")
    if "selectcase(trim(dg_fragment_w90_initial_projection))case('scdm','spectral','random')" not in compact:
        problems.append("fragment W90 projection validation was not preserved")
    return problems


assert not violations(SOURCE), violations(SOURCE)

mutated = SOURCE.replace("case('spectral','random')", "case('spectral','random','invalid')", 1)
assert any("spectral|random" in item for item in violations(mutated))

print("PASS global and fragment W90 initial projection validation with mutation fixture")
