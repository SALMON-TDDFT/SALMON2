#!/usr/bin/env python3
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
FILES = [
    "src/common/dg_hybrid_continuation_residuals.f90",
    "src/common/dg_hybrid_sparse_operators.f90",
    "src/gs/dc/dg_hybrid_continuation_controller.f90",
    "src/gs/dc/dg_hybrid_continuation_scf.f90",
    "src/gs/dc/dg_hybrid_continuation_state.f90",
    "src/gs/dc/dg_hybrid_sipg_operator.f90",
    "tests/dg/test_dg_hybrid_continuation_controller_mpi.f90",
    "tests/dg/test_dg_hybrid_continuation_residuals_mpi.f90",
    "tests/dg/test_dg_hybrid_continuation_scf_mpi.f90",
    "tests/dg/test_dg_hybrid_continuation_state_mpi.f90",
    "tests/dg/test_dg_hybrid_sipg_operator_mpi.f90",
]

failures = []
for relative in FILES:
    path = ROOT / relative
    for number, line in enumerate(path.read_text().splitlines(), 1):
        lowered = line.lower()
        if len(line) > 132:
            failures.append(f"{relative}:{number}: line has {len(line)} columns")
        if re.match(r"^\s*use\s+mpi\s*$", lowered):
            failures.append(f"{relative}:{number}: unrestricted MPI import")
        if re.search(r"::[^!]*\b(comm|rank)\b", lowered):
            failures.append(f"{relative}:{number}: nonconforming parallel identifier")

assert not failures, "\n" + "\n".join(failures)
print("PASS DG continuation coding rules")
