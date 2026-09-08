#!/usr/bin/env python3
from pathlib import Path

source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()
start = source.find("subroutine write_wannier_amn_file_pseudo_channels")
end = source.find("end subroutine write_wannier_amn_file_pseudo_channels", start)
body = source[start:end]

required = [
    "call sawf_projection_shell_lmax(dc%system_tot%nion, nproj",
    "if(complete_shell) then",
    "if(proj_l_raw(ip_raw) > projection_lmax) cycle",
    "selected_count",
    "pseudo_channels selection: complete atom-major shell",
    "call select_top_projectors",
]
missing = [token for token in required if token not in body]
if missing:
    raise SystemExit("complete pseudo-channel shell selection missing: " + ", ".join(missing))

print("PASS complete s+p or s+p+d pseudo-channel shells override projectability ranking")
