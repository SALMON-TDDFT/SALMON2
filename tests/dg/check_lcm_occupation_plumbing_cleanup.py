#!/usr/bin/env python3
"""Reject obsolete occupation-weight plumbing in LCM modules."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
FILES = (
    "src/rt/rt_local_chern_marker.f90",
    "src/rt/rt_local_chern_marker_soi.f90",
)


def main() -> None:
    for relative_path in FILES:
        source = (ROOT / relative_path).read_text()
        for dead_symbol in (
            "local_occ_w",
            "local_occ_weight",
            "local_occ_index",
            "local_occ_global_io",
        ):
            if re.search(rf"\b{dead_symbol}\b", source, re.IGNORECASE):
                raise AssertionError(
                    f"{relative_path}: obsolete symbol remains: {dead_symbol}"
                )

        cache = re.search(
            r"subroutine\s+build_occ_distribution_cache\b(?P<body>.*?)"
            r"end\s+subroutine\s+build_occ_distribution_cache",
            source,
            re.IGNORECASE | re.DOTALL,
        )
        if cache is None:
            raise AssertionError(f"{relative_path}: distribution cache not found")
        if re.search(r"\bocc_w0\b|\blocal_w0\b", cache.group("body"), re.IGNORECASE):
            raise AssertionError(
                f"{relative_path}: distribution cache still transports weights"
            )
        if re.search(
            r"validate_sharp_occupation_contract\s*\([^\n]*\bocc_w\b",
            source,
            re.IGNORECASE,
        ) is None:
            raise AssertionError(
                f"{relative_path}: sharp validator no longer receives occ_w"
            )

    print("LCM occupation plumbing cleanup: PASS")


if __name__ == "__main__":
    main()
