#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
cmake = (ROOT / "src/rt/CMakeLists.txt").read_text().lower()
route = (ROOT / "tests/dg/check_dg_wpw_scf_route.py").read_text().lower()

assert "rt_dg_wpw_distributed_layout.f90" not in cmake
assert not (ROOT / "src/rt/dg/rt_dg_wpw_distributed_layout.f90").exists()
assert "rt_dg_wpw_column_layout.f90" in cmake
assert "rt_dg_wpw_column_layout.f90" in route
print("PASS production route excludes transitional G-only distributed layout")
