#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
src = (root / "src/rt/dg/rt_dg_wpw_length_gauge.f90").read_text()
cmake = (root / "src/rt/CMakeLists.txt").read_text()
required = ["validate_wpw_position_operator", "build_wpw_length_gauge_hamiltonian",
            "evaluate_wpw_polarization", "update_wpw_polarization_branch",
            "h=h0+field*position", "jz=(pz-previous_pz)/dt"]
for token in required:
    assert token in src, token
assert "rt_dg_wpw_length_gauge.f90" in cmake
for forbidden in ("pp_face", "interface_correction", "get_environment_variable"):
    assert forbidden not in src.lower(), forbidden
print("WPW length-gauge source contract: PASS")
