from pathlib import Path


root = Path(__file__).resolve().parents[2]
global_source = (root / "src/io/salmon_global.f90").read_text().lower()
input_source = (root / "src/io/inputoutput.f90").read_text().lower()
lcfo = (root / "src/gs/dc/lcfo_flux.f90").read_text().lower()

name = "wannier_sawf_initial_wavefunction_directory"
assert name in global_source
assert name in input_source
assert "call comm_bcast(wannier_sawf_initial_wavefunction_directory" in input_source

assert "conventional initial wavefunction seed is deferred until pre-diagonalization support" in input_source

writer = lcfo.split("subroutine write_wannier_seed_files", 1)[1].split(
    "end subroutine write_wannier_seed_files", 1
)[0]
assert "call prepare_sawf_initial_wavefunction_seed" not in writer
assert "call apply_sawf_initial_wavefunction_seed_to_lcfo" not in writer

dmn = lcfo.split("subroutine generate_sawf_dmn", 1)[1].split(
    "end subroutine generate_sawf_dmn", 1
)[0]
assert "load_sawf_conventional_reference_states" not in dmn
assert "build_sawf_dmn_operation_global_states" not in dmn
assert "prepare_sawf_fragment_state_cache" in dmn
assert "subroutine prepare_sawf_initial_wavefunction_seed" not in lcfo
assert "subroutine apply_sawf_initial_wavefunction_seed_to_lcfo" not in lcfo
assert "subroutine load_sawf_conventional_reference_states" not in lcfo
assert "subroutine validate_sawf_conventional_reference_receipt" not in lcfo

print("PASS conventional seed fails closed until a coherent pre-diagonalization implementation exists")
