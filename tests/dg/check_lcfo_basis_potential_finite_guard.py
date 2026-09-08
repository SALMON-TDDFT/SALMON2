from pathlib import Path


source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()
assert "call check_lcfo_basis_potential_inputs_finite" in source
assert "[DC-LCFO-INPUT-LOCAL] label=basis" in source
assert "[DC-LCFO-INPUT-LOCAL] label=vlocal" in source
main = source.split('write(*,*) "start DC-LCFO-Flux"', 1)[1]
assert main.index("call check_lcfo_basis_potential_inputs_finite") < main.index("call hpsi_basis")
print("PASS LCFO basis and local potential are diagnosed before Hamiltonian construction")
