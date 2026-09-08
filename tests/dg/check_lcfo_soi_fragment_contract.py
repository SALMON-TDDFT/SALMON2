#!/usr/bin/env python3
"""Protect optimized extents, topology, and same-dvec tags in LCFO-SOI."""

from pathlib import Path
import re


source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_soi.f90").read_text().lower()


def body(name: str) -> str:
    match = re.search(rf"subroutine {name}\b(.*?)end subroutine {name}", source, re.S)
    assert match, f"missing {name}"
    return match.group(1)


init = body("init_lcfo")
basis = body("calc_basis")
hamiltonian = body("calc_hamiltonian_matrix")
output = body("output")

assert "find_fragment_face_neighbor" in source
assert "fragment_direction_tag" in source
assert "call find_fragment_face_neighbor" in init
assert "fragment_extent(:,ifrag)" in init
assert "matches_dst/=1.or.matches_src/=1" in init
assert "dc%nxyz_domain" not in basis
assert "call get_fragment_domain(dc, dc%i_frag, nxyz_domain)" in basis
assert "l = nxyz_domain" in hamiltonian
assert "dc%nxyz_domain" not in hamiltonian
assert re.search(r"fragment_direction_tag\s*\(\s*dc%i_frag\s*,\s*dc%n_frag\s*,\s*halo\(i_halo\)%dvec", hamiltonian)
assert re.search(r"fragment_direction_tag\s*\(\s*halo\(i_halo\)%ifrag_src\s*,\s*dc%n_frag\s*,\s*halo\(i_halo\)%dvec", hamiltonian)
assert "call get_fragment_domain(dc, dc%i_frag, nxyz_domain)" in output

print("PASS LCFO-SOI optimized extents, topology, and same-dvec tags")
