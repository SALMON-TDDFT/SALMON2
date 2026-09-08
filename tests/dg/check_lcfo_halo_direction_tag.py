#!/usr/bin/env python3
"""Protect same-dvec LCFO basis-halo tags and their portable bound."""

from pathlib import Path
import re


root = Path(__file__).resolve().parents[2]
lcfo = (root / "src/gs/dc/lcfo.f90").read_text().lower()
geometry = (root / "src/gs/dc/dc_fragment_geometry.f90").read_text().lower()
match = re.search(r"subroutine calc_hamiltonian_matrix\b(.*?)end subroutine calc_hamiltonian_matrix", lcfo, re.S)
assert match, "missing LCFO Hamiltonian communication"
body = match.group(1)
assert re.search(r"fragment_direction_tag\s*\(\s*dc%i_frag\s*,\s*dc%n_frag\s*,\s*halo\(i_halo\)%dvec", body)
assert re.search(r"fragment_direction_tag\s*\(\s*halo\(i_halo\)%ifrag_src\s*,\s*dc%n_frag\s*,\s*halo\(i_halo\)%dvec", body)
assert "itag_send = dc%i_frag" not in body
assert "itag_recv = halo(i_halo)%ifrag_src" not in body
assert "upper_bound = 32767" in geometry
assert "candidate_tag = int(direction_code,8)*int(n_frag,8)" in geometry
assert "candidate_tag > int(upper_bound,8)" in geometry

print("PASS LCFO sender+same-direction MPI tag contract")
