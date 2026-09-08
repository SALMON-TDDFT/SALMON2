#!/usr/bin/env python3
"""Protect upstream BLAS dot products while retaining DG GS control flow."""

from pathlib import Path
import re


source = (Path(__file__).resolve().parents[2] / "src/gs/gram_schmidt_orth.f90").read_text()
upper = source.upper()

assert upper.count("ZDOTC") == 11, "expected the complete upstream ZDOTC branches"
assert upper.count("COMPLEX(8) :: ZDOTC") == 2, "each complex BLAS routine needs an explicit ZDOTC type"
assert "SUM(CONJG(WF_BLOCK" not in upper, "complex small-block paths must use upstream BLAS ZDOTC"
assert re.search(r"subroutine gram_schmidt_col_cblas\b.*?complex\(8\) :: ZDOTC", source, re.I | re.S)
assert re.search(r"subroutine gram_schmidt_col_cblas_old\b.*?complex\(8\) :: ZDOTC", source, re.I | re.S)
for dg_token in (
    "wf_block_send",
    "if(idiv3 == 0)",
    "call comm_summation(umat_tmp, umat",
):
    assert dg_token.lower() in source.lower(), f"DG Gram-Schmidt behavior lost: {dg_token}"

print("PASS upstream ZDOTC operations and DG Gram-Schmidt flow")
