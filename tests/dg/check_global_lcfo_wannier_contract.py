#!/usr/bin/env python3
"""Source contract for the coherent global LCFO Wannier input space."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
LCFO = (ROOT / "src/gs/dc/lcfo.f90").read_text()
MAIN = (ROOT / "src/gs/main_dft.f90").read_text()


def require(pattern: str, text: str, message: str) -> None:
    if not re.search(pattern, text, re.IGNORECASE | re.DOTALL):
        raise SystemExit(f"FAIL global LCFO Wannier contract: {message}")


require(
    r"subroutine\s+dc_lcfo\s*\(.*?retained_count.*?retained_box_contribution"
    r".*?retained_occupations",
    LCFO,
    "dc_lcfo must expose retained rank, reconstructed columns, and occupations",
)
require(
    r"integer\s*,\s*intent\(in\)\s*,\s*optional\s*::\s*retained_count",
    LCFO,
    "retained_count must be an explicit optional input",
)
require(
    r"coefficient_count\s*=\s*max\s*\(\s*dc%nstate_tot\s*,\s*retained_count\s*\)",
    LCFO,
    "EigenExa/LAPACK must retain every requested coefficient column",
)
require(
    r"present\s*\(\s*retained_count\s*\).*?present\s*\(\s*retained_box_contribution\s*\)"
    r".*?present\s*\(\s*retained_occupations\s*\)",
    LCFO,
    "all retained outputs must be requested as one indivisible contract",
)
require(
    r"ieee_is_finite\s*\(.*?coef_wf",
    LCFO,
    "retained LCFO coefficients must have a finite-value gate",
)
require(
    r"call\s+dc_lcfo\s*\(.*?retained_count\s*=\s*ntarget.*?"
    r"retained_box_contribution\s*=.*?retained_occupations\s*=",
    MAIN,
    "the production adapter must request one coherent retained LCFO space",
)
require(
    r"size\s*\(\s*retained_occupations\s*\)\s*/=\s*retained_count",
    LCFO,
    "occupation output extent must be checked against retained rank",
)
require(
    r"esp_tot\s*\(\s*istate\s*,\s*1\s*\).*?temperature",
    LCFO,
    "fractional occupations must be derived from the retained LCFO eigenvalues",
)
require(
    r"sum\s*\(\s*retained_occupations\s*\).*?electron_count",
    LCFO,
    "retained occupations must conserve the full-system electron count",
)
if re.search(r"retained_occupations.*?system%rocc|system%rocc.*?retained_occupations", LCFO, re.I | re.S):
    raise SystemExit(
        "FAIL global LCFO Wannier contract: fragment-orbital occupations were copied onto LCFO eigenstates"
    )

for obsolete in ("occupied_count", "occupied_box_contribution"):
    if re.search(rf"subroutine\s+dc_lcfo\s*\(.*?\b{obsolete}\b", LCFO, re.I | re.S):
        raise SystemExit(
            "FAIL global LCFO Wannier contract: occupied-only dc_lcfo API remains: "
            + obsolete
        )

print("PASS coherent global LCFO Wannier input contract")
