#!/usr/bin/env python3
"""Check the reduced-memory Si8/MPI8 overlapping-Wannier fixture."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/dg/data/si8_overlapping_wannier"


def require(pattern: str, text: str, description: str) -> None:
    if not re.search(pattern, text, flags=re.IGNORECASE | re.MULTILINE):
        raise RuntimeError(f"Si8 fixture lacks {description}")


def main() -> int:
    input_path = FIXTURE / "inputfile.in"
    atom_path = FIXTURE / "atom.dat"
    if not input_path.is_file() or not atom_path.is_file():
        raise RuntimeError("Si8 fixture files are missing")
    text = input_path.read_text()
    atoms = [line for line in atom_path.read_text().splitlines() if line.strip()]
    if re.search(r"@[A-Z_]+@", text):
        raise RuntimeError("Si8 fixture has unresolved placeholders")
    require(r"\bizatom\s*\(\s*1\s*\)\s*=\s*14\b", text, "silicon atomic number")
    require(r"\bnatom\s*=\s*8\b", text, "eight atoms")
    require(r"\bnelec\s*=\s*32\b", text, "32 electrons")
    require(r"\bnum_rgrid\s*\(\s*1\s*:\s*3\s*\)\s*=\s*16\s*,\s*16\s*,\s*16", text, "16-cubed grid")
    require(r"\bnum_fragment\s*\(\s*1\s*:\s*3\s*\)\s*=\s*2\s*,\s*2\s*,\s*2", text, "2x2x2 fragments")
    require(r"\bnproc_rgrid_tot\s*\(\s*1\s*:\s*3\s*\)\s*=\s*2\s*,\s*2\s*,\s*2", text, "2x2x2 rank grid")
    require(r"\byn_dg_dc_overlapping_wannier\s*=\s*['\"]y['\"]", text, "overlapping-Wannier route")
    require(r"\bfile_pseudo\s*\(\s*1\s*\)\s*=\s*['\"]Si_rps\.dat['\"]", text, "Si pseudopotential")
    if len(atoms) != 8 or any(not line.lstrip().startswith("'Si'") for line in atoms):
        raise RuntimeError("Si8 atom.dat must contain exactly eight silicon atoms")
    print("Si8 overlapping-Wannier fixture contract: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
