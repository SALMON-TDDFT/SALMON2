from pathlib import Path


source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()
body = source.split("subroutine generate_sawf_dmn", 1)[1].split(
    "end subroutine generate_sawf_dmn", 1
)[0]

assert "SAWF split-fragment symmetry requires a constructed global symmetry-closed basis" in body
guard = body.index("SAWF split-fragment symmetry requires a constructed global symmetry-closed basis")
local_build = body.index("call build_sawf_dmn_operation_fragment_local")
assert guard < local_build
assert "build D_band from split fragment blocks" not in body

print("PASS split-fragment symmetry cannot masquerade as a global closed SAWF basis")
