#!/usr/bin/env python3
"""Regression gate for the typed interface-continuation diagnostic record."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/main_dft.f90").read_text().lower()
routine = source.split(
    "subroutine record_dg_hybrid_interface_continuation_diagnostic", 1
)[1].split("end subroutine record_dg_hybrid_interface_continuation_diagnostic", 1)[0]
record = re.search(
    r"write\s*\(\s*\*\s*,\s*'([^']+)'\s*\).*?"
    r"\[dg-hybrid-continuation\]\s+lambda=",
    routine,
    re.S,
)
assert record, "missing interface-continuation diagnostic write"
assert record.group(1).startswith("(2(a,es24.16)"), (
    "diagnostic format has an extra leading character descriptor"
)

print("interface-continuation diagnostic format: PASS")
