#!/usr/bin/env python3
"""Require LCFO core-array sites to use optimized per-fragment extents."""

from pathlib import Path
import re


source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo.f90").read_text()
geometry = (Path(__file__).resolve().parents[2] / "src/gs/dc/dc_fragment_geometry.f90").read_text().lower()


def body(name: str) -> str:
    match = re.search(rf"subroutine {name}\b(.*?)end subroutine {name}", source, re.I | re.S)
    assert match, f"missing {name}"
    return match.group(1).lower()


calc_basis = body("calc_basis")
init_lcfo = body("init_lcfo")
hamiltonian = body("calc_hamiltonian_matrix")
output = body("output")
test_write = body("test_write_psi")

assert "dc%nxyz_domain" not in calc_basis, "calc_basis still indexes optimized storage by nominal extent"
assert "call get_fragment_domain(dc, dc%i_frag, nxyz_domain)" in calc_basis
assert "l = nxyz_domain" in hamiltonian, "Hamiltonian block still truncates/overruns optimized basis"
assert "dc%nxyz_domain" not in hamiltonian
assert output.index("call get_fragment_domain(dc, dc%i_frag, nxyz_domain)") < output.index("if(dc%id_frag==0)"), (
    "every rank must initialize the output extent before buffered indexing"
)
assert "dc%nxyz_domain" not in test_write, "diagnostic reconstruction still uses nominal extent"
assert re.search(r"call\s+get_fragment_domain\s*\(\s*dc\s*,\s*dc%i_frag\s*,\s*nxyz_domain\s*\)", test_write)
assert "call find_fragment_face_neighbor" in init_lcfo
assert "*nxyz_domain" not in init_lcfo, "neighbor origins must not assume the current fragment width"
assert "fragment_extent(:,ifrag)" in init_lcfo
assert "matches_dst/=1.or.matches_src/=1" in init_lcfo
assert "subroutine find_fragment_face_neighbor" in geometry
assert "fragment_extent(axis,candidate)" in geometry, "negative-face matching must use candidate extent"
assert "modulo(" in geometry, "periodic face matching must use nonnegative modulo"

print("PASS LCFO uses per-fragment extents for core arrays and periodic face topology")
