#!/usr/bin/env python3
"""Guard collective failure for ill-conditioned LCM dual overlaps."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]


def require(source: str, pattern: str, message: str) -> None:
    if re.search(pattern, source, re.IGNORECASE | re.DOTALL) is None:
        raise AssertionError(message)


def check_source(relative_path: str) -> None:
    source = (ROOT / relative_path).read_text()
    require(
        source,
        r"subroutine\s+build_transposed_inverse_coefficients_checked\s*\("
        r"\s*a\s*,\s*label\s*,\s*inversion_ok\s*,\s*rcond_out\s*\)",
        f"{relative_path}: root inverse helper does not return status and rcond",
    )
    require(
        source,
        r"if\s*\(\s*rcond\s*<\s*rcond_warn\s*\)\s*then.*?"
        r"inversion_ok\s*=\s*\.false\..*?return",
        f"{relative_path}: low-rcond path does not fail before zgetrs",
    )
    require(
        source,
        r"call\s+comm_bcast\s*\(\s*inversion_ok\s*,\s*info%icomm_o\s*,\s*0\s*\)",
        f"{relative_path}: inversion status is not broadcast on icomm_o",
    )
    require(
        source,
        r"call\s+comm_bcast\s*\(\s*rcond_root\s*,\s*info%icomm_o\s*,\s*0\s*\)",
        f"{relative_path}: reciprocal condition estimate is not broadcast",
    )
    require(
        source,
        r"if\s*\(\s*\.not\.\s*inversion_ok\s*\)\s*then.*?"
        r"stop\s+['\"]LCM dual overlap is ill-conditioned['\"]",
        f"{relative_path}: orbital ranks do not stop collectively",
    )


def main() -> None:
    check_source("src/rt/rt_local_chern_marker.f90")
    check_source("src/rt/rt_local_chern_marker_soi.f90")
    print("LCM collective condition guard: PASS")


if __name__ == "__main__":
    main()
