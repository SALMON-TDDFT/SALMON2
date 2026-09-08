#!/usr/bin/env python3
"""Guard against mixing Wannier90 gauge data with bond-center labels."""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src" / "gs" / "dc" / "lcfo_flux.f90"


def main() -> int:
    text = SOURCE.read_text()
    failures = []

    if re.search(r"write\s*\(\s*iunit\s*\)\s*export_center_bohr", text, re.IGNORECASE):
        failures.append("Wannier90 global basis must not export substituted bond-center coordinates.")

    if re.search(r"export_center_bohr\s*\([^\\n]*\)\s*=\s*owner_center_bohr", text, re.IGNORECASE):
        failures.append("Bond-center coordinates are being assigned to the exported W90 center array.")

    if "exported center source=bond_centers" in text:
        failures.append("Log message still advertises bond-center export for W90 gauge data.")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1

    print("PASS: W90 global basis export keeps W90 centers with W90 gauge data.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
