#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
types = root / "src/rt/dg/rt_dg_nodal_types.f90"
cmake = (root / "src/rt/CMakeLists.txt").read_text()

assert types.exists(), "missing nodal real-space DG types"
body = types.read_text()
assert "type :: nodal_face_info" in body
assert "integer :: axis" in body
assert "integer :: side" in body
assert "integer :: neighbor_fragment" in body
assert "complex(8), allocatable :: psi_core" in body
assert "complex(8), allocatable :: recv_value" in body
assert "complex(8), allocatable :: recv_normal" in body
assert "psi_buffer" not in body, "buffer must be a halo, not an independently propagated DOF"
assert "dg/rt_dg_nodal_types.f90" in cmake
print("PASS nodal DG stores core DOFs and explicit face halos")
