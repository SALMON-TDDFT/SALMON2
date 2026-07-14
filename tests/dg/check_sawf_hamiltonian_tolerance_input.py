from pathlib import Path

root = Path(__file__).resolve().parents[2]
global_source = (root / "src/io/salmon_global.f90").read_text().lower()
input_source = (root / "src/io/inputoutput.f90").read_text().lower()
flux = (root / "src/gs/dc/lcfo_flux.f90").read_text().lower()
dmn = (root / "src/gs/dc/lcfo_wannier_sawf_dmn.f90").read_text().lower()

name = "wannier_sawf_hamiltonian_tolerance"
assert name in global_source
assert name in input_source
assert f"call comm_bcast({name}" in input_source
assert name in flux
assert "hamiltonian_tolerance=" in flux
assert "writer%hamiltonian_tolerance" in dmn
assert "hres>hamiltonian_tolerance_effective" in dmn
generator = flux.split("subroutine generate_sawf_dmn", 1)[1].split(
    "end subroutine generate_sawf_dmn", 1
)[0]
assert "hamiltonian_tolerance=" in generator
assert "writer%hamiltonian_tolerance" in generator

print("PASS SAWF Hamiltonian covariance has an explicit namelist tolerance")
