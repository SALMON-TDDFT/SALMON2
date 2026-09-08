#!/usr/bin/env python3
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text()

generator = source.split("subroutine generate_sawf_dmn", 1)[1].split(
    "end subroutine generate_sawf_dmn", 1
)[0]
required = [
    "symmetry_fragment_maps",
    "call build_sawf_closed_fragment_seed_basis",
    "type(t_sawf_closed_basis)",
]
for token in required:
    if token not in generator:
        raise SystemExit(f"SAWF generator is missing {token}")

prepare = generator.index("call prepare_sawf_fragment_state_cache")
closed = generator.index("call build_sawf_closed_fragment_seed_basis")
d_band = generator.index("call build_sawf_dmn_operation_fragment_local")
if not prepare < closed < d_band:
    raise SystemExit("SAWF closed seed must be built after cache preparation and before D_band")

builder = source.split("subroutine build_sawf_closed_fragment_seed_basis", 1)[1].split(
    "end subroutine build_sawf_closed_fragment_seed_basis", 1
)[0]
for token in [
    "build_sawf_fragment_buffer_point_map",
    "append_sawf_mapped_basis",
    "build_sawf_closed_core_buffer_basis",
    "findloc(fragment_maps(:,iop)==dc%i_frag",
    "source_entry%buffer_basis",
]:
    if token not in builder.replace(" ", "") if "findloc" in token else token not in builder:
        raise SystemExit(f"SAWF closed fragment builder is missing {token}")

print("PASS SAWF closed seed collects core and buffer symmetry orbits")
