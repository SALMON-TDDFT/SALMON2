#!/usr/bin/env python3
"""Contracts for upstream DC input correctness backports."""

import re
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = (root / "src/io/inputoutput.f90").read_text()

# The integrated input validator contains earlier nested yn_dc checks for DG.
# Match the actual restart-rewrite clause rather than truncating at the first
# unrelated nested ``end if``.
assert re.search(r"if\s*\(\s*write_gs_restart_data\s*/=\s*'no'\s*\)\s*then.*?"
                 r"write_gs_restart_data\s*=\s*'no'", source, re.I | re.S), (
    "DC must disable the incompatible conventional GS restart writer"
)
print("upstream DC input backports: PASS")
