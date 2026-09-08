from pathlib import Path


root = Path(__file__).resolve().parents[2]
global_text = (root / "src/io/salmon_global.f90").read_text().lower()
input_text = (root / "src/io/inputoutput.f90").read_text().lower()
restart_text = (root / "src/io/checkpoint_restart.f90").read_text().lower()

name = "yn_reset_occupation_restart"
assert name in global_text
assert name in input_text
assert f"call comm_bcast({name}" in input_text
assert name in restart_text
assert "flag_read_occ  = .false." in restart_text.split(name, 1)[1]

print("PASS restart occupation reset is an explicit namelist-controlled operation")
