#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
halo = root / "src/rt/dg/rt_dg_nodal_halo.f90"
cmake = (root / "src/rt/CMakeLists.txt").read_text()

assert halo.exists(), "missing nodal face-halo topology"
body = halo.read_text()
assert "do axis = 1, 3" in body
assert "do side = -1, 1, 2" in body
assert "iface = nodal_face_slot(axis, side)" in body
assert "faces(iface)%axis = axis" in body
assert "faces(iface)%side = side" in body
assert "neighbor_fragment_from_coords" in body
assert "modulo(" in body, "periodic fragment topology must wrap"
assert "ifrag = 1 + coords(3) + num_fragment(3)" in body
assert "coords(2) + num_fragment(2) * coords(1)" in body
assert "if (faces(iface)%axis /= axis .or. faces(iface)%side /= side)" in body
assert "dg/rt_dg_nodal_halo.f90" in cmake
print("PASS nodal halo keeps +/- periodic faces distinct")
