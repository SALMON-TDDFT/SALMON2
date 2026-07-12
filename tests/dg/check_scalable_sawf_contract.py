from pathlib import Path

root = Path(__file__).resolve().parents[2]
contract = (root / "docs/plans/2026-07-12-scalable-sawf-contract.md").read_text().lower()
global_input = (root / "src/io/salmon_global.f90").read_text().lower()
inputoutput = (root / "src/io/inputoutput.f90").read_text().lower()

required = [
    "monolithic global sawf", "actual supercell symmetry group", "parent-crystal",
    "d_band", "d_wann", "complete projection shell", "representative environment",
    "defect-local regeneration", "buffer convergence", "unitary procrustes",
    "template provenance", "cache invalidation", "dg face block",
]
for token in required:
    assert token in contract, f"contract is missing {token}"
assert "site_symmetry=true" in contract and "not an acceptance" in contract

controls = [
    "wannier_sawf_generation", "wannier_sawf_symmetry_scope",
    "wannier_sawf_parent_symmetry_file", "wannier_sawf_buffer_steps",
    "wannier_sawf_gauge_tolerance", "wannier_sawf_buffer_tolerance",
    "wannier_sawf_equivalence_tolerance", "wannier_sawf_cache_directory",
]
for name in controls:
    assert name in global_input, f"missing namelist storage for {name}"
    assert name in inputoutput, f"missing namelist route for {name}"
assert "get_environment_variable" not in inputoutput.lower()
print("PASS scalable SAWF contract and namelist route")
