#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = (root / "src/gs/dc/lcfo_flux.f90").read_text().lower()
body = source.split("subroutine prepare_sawf_closed_seed_eigensystem", 1)[1]
body = body.split("end subroutine prepare_sawf_closed_seed_eigensystem", 1)[0]
load = body.split("call load_sawf_symmetry_file", 1)[1]
before_reduce = load.split("call comm_get_max(failure", 1)[0]

assert "local_ok" in before_reduce
assert "message" in before_reduce
assert "symmetry_filename" in before_reduce
assert "dc%id_tot" in before_reduce
assert "sawf closed-seed symmetry load" in before_reduce

print("PASS SAWF symmetry-load failure is diagnosed before collective reduction")
