#!/usr/bin/env python3
"""Guard the sharp all-electron Local Chern Marker contract."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
SCALAR = (ROOT / "src/rt/rt_local_chern_marker.f90").read_text()
SOI = (ROOT / "src/rt/rt_local_chern_marker_soi.f90").read_text()


def require(source: str, pattern: str, message: str) -> None:
    if re.search(pattern, source, re.IGNORECASE | re.DOTALL) is None:
        raise AssertionError(message)


def reject(source: str, pattern: str, message: str) -> None:
    if re.search(pattern, source, re.IGNORECASE | re.DOTALL) is not None:
        raise AssertionError(message)


def main() -> None:
    require(
        SCALAR,
        r"if\s*\(\s*system%nspin\s*==\s*1\s*\)\s*then\s*"
        r"electron_multiplicity\s*=\s*2\.0d0",
        "scalar LCM does not assign multiplicity 2 for nspin=1",
    )
    require(
        SCALAR,
        r"else\s*\n\s*electron_multiplicity\s*=\s*1\.0d0",
        "scalar LCM does not assign multiplicity 1 for explicit spin channels",
    )
    require(
        SOI,
        r"electron_multiplicity\s*=\s*1\.0d0",
        "SOI LCM does not assign spinor multiplicity 1",
    )
    require(
        SCALAR,
        r"call\s+validate_sharp_occupation_contract\s*\([^\n]*"
        r"electron_multiplicity",
        "scalar LCM does not validate sharp occupations",
    )
    require(
        SOI,
        r"call\s+validate_sharp_occupation_contract\s*\([^\n]*"
        r"electron_multiplicity",
        "SOI LCM does not validate sharp occupations",
    )
    require(
        SCALAR,
        r"marker_local\s*\([^\n]*&\s*\n[^\n]*electron_multiplicity",
        "scalar marker contribution does not include electron multiplicity",
    )
    require(
        SOI,
        r"marker_local\s*\([^\n]*&\s*\n[^\n]*electron_multiplicity",
        "SOI marker contribution does not include electron multiplicity",
    )
    reject(
        SCALAR,
        r"sqrt\s*\(\s*max\s*\(\s*0\.0d0\s*,\s*local_occ_w",
        "scalar LCM still applies occupation weights before Lowdin",
    )
    reject(
        SOI,
        r"sqrt\s*\(\s*max\s*\(\s*0\.0d0\s*,\s*local_occ_w",
        "SOI LCM still applies occupation weights before Lowdin",
    )
    require(
        SCALAR + SOI,
        r"unsupported\s+fractional\s+occupation",
        "LCM has no fail-closed diagnostic for fractional occupation",
    )
    print("LCM all-electron occupation contract: PASS")


if __name__ == "__main__":
    main()
