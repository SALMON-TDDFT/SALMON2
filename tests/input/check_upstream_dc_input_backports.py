#!/usr/bin/env python3
"""Contracts for upstream DC input correctness backports."""

import re
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = (root / "src/io/inputoutput.f90").read_text()

dc = re.search(r"if\s*\(\s*yn_dc\s*==\s*'y'\s*\)\s*then(?P<body>.*?)end\s+if",
               source, re.I | re.S)
assert dc, "DC input validation block is missing"
body = dc.group("body")
assert re.search(r"write_gs_restart_data\s*/=\s*'no'.*?"
                 r"write_gs_restart_data\s*=\s*'no'", body, re.I | re.S), (
    "DC must disable the incompatible conventional GS restart writer"
)
print("upstream DC input backports: PASS")
