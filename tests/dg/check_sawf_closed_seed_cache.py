#!/usr/bin/env python3
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text()

entry = source.split("type :: t_sawf_fragment_state_entry", 1)[1].split(
    "end type t_sawf_fragment_state_entry", 1
)[0]
if not re.search(r"allocatable\s*::[^\n]*buffer_basis\s*\(", entry, re.I):
    raise SystemExit("SAWF fragment cache does not retain the buffered basis")
if not re.search(r"integer\s*::[^\n]*buffer_shape", entry, re.I):
    raise SystemExit("SAWF fragment cache does not retain buffer geometry")

loader = source.split("subroutine load_sawf_fragment_state_entry", 1)[1].split(
    "end subroutine load_sawf_fragment_state_entry", 1
)[0]
required = [
    "read_fragment_lcfo_buffer_seed_for_wannier_import",
    "nxyz_buffer_read",
    "buffer_basis",
    "n_basis_buffer",
]
for token in required:
    if token not in loader:
        raise SystemExit(f"SAWF fragment cache loader is missing {token}")
if "any(nxyz_buffer_read/=dc%nxyz_buffer)" not in loader.replace(" ", ""):
    raise SystemExit("SAWF fragment cache does not validate current buffer width")
if "n_basis_buffer/=n_basis_frag" not in loader.replace(" ", ""):
    raise SystemExit("SAWF fragment cache does not bind core and buffer basis counts")

clearer = source.split("subroutine clear_sawf_fragment_state_entry", 1)[1].split(
    "end subroutine clear_sawf_fragment_state_entry", 1
)[0]
if "deallocate(entry%buffer_basis)" not in clearer.replace(" ", ""):
    raise SystemExit("SAWF fragment cache does not release the buffered basis")

print("PASS SAWF current-run cache binds core and buffered fragment bases")
