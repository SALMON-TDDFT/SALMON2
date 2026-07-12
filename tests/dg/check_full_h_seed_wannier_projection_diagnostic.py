from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "rt" / "dg" / "rt_dg_fragment.f90"


text = SRC.read_text()

required = [
    "[DG-FULL-H-SEED-XI-PROJ]",
    "proj_norm_min_occ",
    "proj_defect_max_occ",
    "proj_norm_min_all",
    "proj_defect_max_all",
]

missing = [token for token in required if token not in text]
if missing:
    raise SystemExit("missing Full-H/Wannier projection diagnostic tokens: " + ", ".join(missing))
