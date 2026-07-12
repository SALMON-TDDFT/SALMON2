from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/dc/lcfo_wannier_sawf_orchestrator.f90").read_text().lower()
flux = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text().lower()

required = [
    "type, public :: t_sawf_environment_receipt",
    "build_sawf_environment_execution_plan",
    "validate_sawf_environment_receipts",
    "representative_fragment",
    "operation_index",
    "generated_independently",
    "same_supercell_fingerprint",
]
for token in required:
    assert token in source, f"missing hierarchical orchestrator contract: {token}"

assert "call build_sawf_environment_execution_plan" in flux

print("PASS production hierarchical SAWF execution-plan contract")
