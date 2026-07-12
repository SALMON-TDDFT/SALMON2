#!/usr/bin/env python3
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text()

loader = source.split("subroutine load_sawf_fragment_state_entry", 1)[1].split(
    "end subroutine load_sawf_fragment_state_entry", 1
)[0]
if "sawf_explicit_basis_active" not in loader or "read_sawf_closed_fragment_state_entry" not in loader:
    raise SystemExit("SAWF D_band cache does not switch to the closed seed")

reader = source.split("subroutine read_sawf_closed_fragment_state_entry", 1)[1].split(
    "end subroutine read_sawf_closed_fragment_state_entry", 1
)[0]
for token in [
    "sawf_closed_seed_magic",
    "sawf_closed_seed_version",
    "current_sawf_seed_token",
    "entry%buffer_basis",
    "entry%states=matmul",
]:
    if token.replace(" ", "") not in reader.replace(" ", ""):
        raise SystemExit(f"closed SAWF seed reader is missing {token}")

writer = source.split("subroutine write_sawf_closed_seed_file", 1)[1].split(
    "end subroutine write_sawf_closed_seed_file", 1
)[0]
for token in ["current_sawf_seed_token", "f_basis", "sawf_explicit_buffer", "coef_wf", "esp_tot"]:
    if token not in writer:
        raise SystemExit(f"closed SAWF seed writer is missing {token}")

print("PASS SAWF D_band cache is bound to the closed eigenseed")
