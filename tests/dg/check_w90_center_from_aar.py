from pathlib import Path


src = Path("src/gs/dc/lcfo_flux.f90").read_text()

if "set_wannier_centers_from_position_diagonal" not in src:
    raise SystemExit("missing AA_R-diagonal Wannier center synchronization helper")

if src.count("center source=AA_R diagonal") < 2:
    raise SystemExit("both W90 export and import-only paths must log AA_R center source")

for anchor in [
    "call read_wannier90_global_rmn_gamma_block_import",
    "call read_wannier90_global_rmn_gamma_block(num_wann_chk",
]:
    i = src.find(anchor)
    if i < 0:
        raise SystemExit(f"missing AA_R read anchor: {anchor}")
    j = src.find("set_wannier_centers_from_position_diagonal", i)
    k = src.find("write(iunit) center_bohr", i)
    if j < 0 or k < 0 or j > k:
        raise SystemExit(
            "Wannier centers must be synchronized from AA_R diagonal before "
            f"writing global basis after {anchor}"
        )

print("W90 global-basis centers are synchronized from AA_R diagonals")
