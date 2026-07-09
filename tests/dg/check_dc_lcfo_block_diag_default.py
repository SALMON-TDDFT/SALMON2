from pathlib import Path
import re


root = Path(".")
inputoutput = (root / "src/io/inputoutput.f90").read_text()
salmon_global = (root / "src/io/salmon_global.f90").read_text()
lcfo_flux = (root / "src/gs/dc/lcfo_flux.f90").read_text()

if "yn_dc_lcfo_block_diag_h" not in salmon_global:
    raise SystemExit("yn_dc_lcfo_block_diag_h must be declared in salmon_global")

required_input_patterns = [
    r"namelist/dc/.*yn_dc_lcfo_block_diag_h",
    r"yn_dc_lcfo_block_diag_h\s*=\s*'n'",
    r"call\s+comm_bcast\(yn_dc_lcfo_block_diag_h",
    r"yn_argument_check\(yn_dc_lcfo_block_diag_h\)",
]
for pattern in required_input_patterns:
    if not re.search(pattern, inputoutput, re.S):
        raise SystemExit(f"missing input handling for yn_dc_lcfo_block_diag_h: {pattern}")

if re.search(r"enabled\s*=\s*\.true\.\s*\n\s*env_value\s*=", lcfo_flux):
    raise SystemExit("block-diag H must not be enabled by default")

if not re.search(r"enabled\s*=\s*\(yn_dc_lcfo_block_diag_h\s*==\s*'y'\)", lcfo_flux):
    raise SystemExit("block-diag H default must come from yn_dc_lcfo_block_diag_h")

print("DC-LCFO block-diag H default is namelist-controlled and disabled")
