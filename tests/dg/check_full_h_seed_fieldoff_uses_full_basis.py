from pathlib import Path


src = Path("src/rt/dg/rt_dg_integrator_expdiag.f90").read_text()

full_seed_branch = src.find("if (yn_dg_full_h_eigen_seed == 'y') then")
legacy_phase_branch = src.find("use_full_h_seed_phase =")

if full_seed_branch < 0:
    raise SystemExit("missing full-DG eigen seed production branch")
if legacy_phase_branch < 0:
    raise SystemExit("missing legacy field-free phase branch marker")

if not full_seed_branch < legacy_phase_branch:
    raise SystemExit(
        "field-free full-DG eigen-seed propagation must use "
        "apply_full_h_seed_eigen_exp before the legacy esp-only phase fallback"
    )

branch_end = src.find("use_full_h_seed_phase =", full_seed_branch)
branch = src[full_seed_branch:branch_end]

required = [
    "dg_frag%has_full_h_seed_eigen",
    "dg_frag%has_full_h_seed_xi",
    "call apply_full_h_seed_eigen_exp(state_first, state_last)",
]
for needle in required:
    if needle not in branch:
        raise SystemExit(f"full-DG eigen seed branch is missing: {needle}")
