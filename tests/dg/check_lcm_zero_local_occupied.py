#!/usr/bin/env python3
"""Guard the LCM copy helpers against dummy orbital columns."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]


def check_source(relative_path: str, rank: int) -> None:
    source = (ROOT / relative_path).read_text()
    errors: list[str] = []

    if not re.search(
        r"call\s+copy_occupied_to_temp\s*\([^\n]*nocc_local\s*,\s*zocc\s*\)",
        source,
        re.IGNORECASE,
    ):
        errors.append("copy call does not pass the actual nocc_local")

    helper_match = re.search(
        r"subroutine\s+copy_occupied_to_temp\b(?P<body>.*?)"
        r"end\s+subroutine\s+copy_occupied_to_temp",
        source,
        re.IGNORECASE | re.DOTALL,
    )
    if helper_match is None:
        errors.append("copy_occupied_to_temp helper not found")
    else:
        body = helper_match.group("body")
        if not re.search(r"nocc_local_in", body, re.IGNORECASE):
            errors.append("helper has no explicit local occupied-count argument")
        if not re.search(
            r"nocc_local0\s*=\s*nocc_local_in\b", body, re.IGNORECASE
        ):
            errors.append("helper does not use the explicit occupied count")
        if re.search(
            rf"nocc_local0\s*=\s*size\s*\(\s*zbuf\s*,\s*{rank}\s*\)",
            body,
            re.IGNORECASE,
        ):
            errors.append("helper still treats the dummy array extent as occupied")

    if errors:
        raise AssertionError(f"{relative_path}: " + "; ".join(errors))


def main() -> None:
    check_source("src/rt/rt_local_chern_marker.f90", rank=4)
    check_source("src/rt/rt_local_chern_marker_soi.f90", rank=5)
    print("LCM zero-local-occupied contract: PASS")


if __name__ == "__main__":
    main()
