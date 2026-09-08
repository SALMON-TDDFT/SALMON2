#!/usr/bin/env python3
"""Require production DG rank/fragment validation in input checking."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/io/inputoutput.f90").read_text(errors="replace").lower()


def violations(source: str) -> list[str]:
    start = source.index("subroutine check_bad_input")
    body = source[start : source.index("end subroutine check_bad_input", start)]
    compact = re.sub(r"\s+|&", "", body)
    route = (
        "yn_dg_dc_overlapping_wannier=='y'.or.yn_dg_hybrid_divided_scf=='y'.or."
        "yn_dg_hybrid_continuation_scf=='y'.or.yn_dg_hybrid_scf=='y'"
    )
    blocks = {
        "positive": (
            f"if(({route}).and.any(num_fragment<=0))"
            'callsawf_input_fatal("productionoverlapping/hybriddgrequirespositivenum_fragment")'
        ),
        "rank": (
            f"if(({route}).and.nproc_size_global/=product(num_fragment))"
            'callsawf_input_fatal("productionoverlapping/hybriddgrequiresonempirankperfragment")'
        ),
        "uniform": (
            f"if(({route}).and.yn_dc_fragment_optimization=='y')"
            'callsawf_input_fatal("productionoverlapping/hybriddgrequiresuniformdcfragments")'
        ),
    }
    problems: list[str] = []
    positions: dict[str, int] = {}
    for name, block in blocks.items():
        positions[name] = compact.find(block)
        if positions[name] < 0:
            problems.append(f"early {name} validation is not coupled to the complete route predicate and diagnostic")
    detailed_position = compact.find("if(dg_dc_handoff_min_iter<1)")
    if detailed_position < 0 or any(position < 0 or position > detailed_position for position in positions.values()):
        problems.append("all rank/fragment validations must precede detailed DG controls")
    elif not positions["positive"] < positions["rank"] < positions["uniform"]:
        problems.append("positive, rank, and uniform validations are out of order")
    return problems


assert not violations(SOURCE), violations(SOURCE)

# Independent semantic mutations must invalidate each complete early block.
no_positive = SOURCE.replace("any(num_fragment<=0)", "all(num_fragment>0)", 1)
assert any("early positive validation" in item for item in violations(no_positive))
no_early_rank = SOURCE.replace("nproc_size_global /= product(num_fragment)", "nproc_size_global == product(num_fragment)", 1)
assert any("early rank validation" in item for item in violations(no_early_rank))
no_uniform = SOURCE.replace("yn_dc_fragment_optimization=='y'", "yn_dc_fragment_optimization/='y'", 1)
assert any("early uniform validation" in item for item in violations(no_uniform))

# A late check in the same input routine, and the separate runtime defenses in
# main_dft, cannot satisfy this early complete-block contract.
late_rank = no_early_rank.replace(
    "end subroutine check_bad_input",
    "if(nproc_size_global /= product(num_fragment)) stop 'late mutation'\n  end subroutine check_bad_input",
    1,
)
assert any("early rank validation" in item for item in violations(late_rank))

print("PASS production DG rank/fragment and uniformity checks occur during early input validation")
