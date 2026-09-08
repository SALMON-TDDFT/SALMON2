from pathlib import Path


root = Path(__file__).resolve().parents[2]
global_text = (root / "src/io/salmon_global.f90").read_text()
input_text = (root / "src/io/inputoutput.f90").read_text()
contract = " ".join(
    (root / "docs/plans/2026-07-12-scalable-sawf-contract.md").read_text().lower().split()
)

assert "wannier_sawf_global_reference_source" in global_text
assert "wannier_sawf_global_reference_source" in input_text
assert "call comm_bcast(wannier_sawf_global_reference_source" in input_text
assert "case('lcfo')" in input_text
assert "accepted SAWF global reference source must be lcfo" in input_text
assert "same actual supercell" in contract
assert "conventional checkpoint is not an accepted sawf reference" in contract

print("PASS monolithic SAWF reference source is an explicit namelist contract")
