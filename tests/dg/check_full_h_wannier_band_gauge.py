#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = (root / "src/rt/dg/rt_dg_fragment.f90").read_text()
assert "build_full_h_seed_band_gauge_position_operator" in source
assert "gauge_orth = matmul(gauge,metric_invhalf)" in source
assert "[DG-FULL-H-BAND-GAUGE-PRINCIPAL]" in source
input_source = (root / "src/io/inputoutput.f90").read_text()
assert "yn_dg_full_h_wannier_band_gauge = 'n'" in input_source
print("PASS Full-H Wannier band-gauge position route is wired")
