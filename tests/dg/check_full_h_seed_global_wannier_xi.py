from pathlib import Path


src = Path("src/rt/dg/rt_dg_fragment.f90").read_text()

required = [
    "build_full_h_seed_global_wannier_position_operator",
    "dg_frag%has_global_wannier_position",
    "dg_frag%global_wannier_position",
    "source=global_wannier_position_AA_R",
]
for token in required:
    if token not in src:
        raise SystemExit(f"missing Full-H seed global-Wannier xi support: {token}")

global_branch = src.find("call build_full_h_seed_global_wannier_position_operator()")
sawtooth_source = src.find("source=centered-cell-sawtooth")
if global_branch < 0 or sawtooth_source < 0 or global_branch > sawtooth_source:
    raise SystemExit(
        "Full-H seed xi must prefer the W90 global position matrix before "
        "falling back to centered-cell sawtooth coordinates"
    )

print("Full-H seed xi prefers the W90 global position matrix when available")
