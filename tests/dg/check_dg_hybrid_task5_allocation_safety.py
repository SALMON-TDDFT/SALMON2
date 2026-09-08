#!/usr/bin/env python3
"""Source contract for Task 5 allocation-failure safety."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
COMPLEMENT = ROOT / "src/common/dg_hybrid_wannier_complement.f90"
STREAM = ROOT / "src/gs/dc/dg_hybrid_fragment_basis_stream.f90"


def subroutine_body(source: str, name: str) -> str:
    match = re.search(
        rf"\bsubroutine\s+{name}\b(?P<body>.*?)"
        rf"\bend\s+subroutine\s+{name}\b",
        source,
        re.IGNORECASE | re.DOTALL,
    )
    if match is None:
        raise AssertionError(f"missing subroutine {name}")
    return match.group("body")


def compact_fortran(source: str) -> str:
    source = re.sub(r"!.*", "", source)
    source = re.sub(r"&\s*\n\s*&?", "", source)
    return re.sub(r"\s+", "", source).lower()


def logical_statements(source: str) -> list[str]:
    statements: list[str] = []
    pending = ""
    for raw_line in source.splitlines():
        line = raw_line.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("&"):
            line = line[1:].lstrip()
        continued = line.endswith("&")
        if continued:
            line = line[:-1]
        pending += line
        if not continued:
            statements.append(re.sub(r"\s+", "", pending).lower())
            pending = ""
    if pending:
        statements.append(re.sub(r"\s+", "", pending).lower())
    return statements


def allocation_statement(body: str, anchor: str) -> str:
    candidates = [
        statement
        for statement in logical_statements(body)
        if statement.startswith("allocate(") and anchor in statement
    ]
    if len(candidates) != 1:
        raise AssertionError(
            f"expected one allocation containing {anchor}, found {len(candidates)}"
        )
    return candidates[0]


def require_checked_local_allocation(body: str) -> None:
    compact = compact_fortran(body)
    statement = allocation_statement(body, "projected_pw_buffer(width,npoint)")
    assert "stat=" in statement, "projected PW allocation must use stat="
    allocation = compact.index(statement)
    payload = compact.index("projected_pw_buffer=raw_pw_buffer", allocation)
    failure_path = compact[allocation + len(statement) : payload]
    assert "allocation_status" in failure_path, (
        "projected PW allocation status must be checked before writing the payload"
    )
    assert "message=" in failure_path and "return" in failure_path, (
        "projected PW allocation failure must return an explanatory failure"
    )


def require_collective_allocation_phase(
    body: str, allocation_anchor: str, first_payload_use: str, label: str
) -> None:
    compact = compact_fortran(body)
    statement = allocation_statement(body, allocation_anchor)
    assert "stat=" in statement, f"{label} allocation must use stat="
    allocation = compact.index(statement)
    payload = compact.index(first_payload_use, allocation + len(statement))
    guard = compact[allocation + len(statement) : payload]
    assert "allocation_status" in guard, (
        f"{label} allocation status must be converted to a rank-local failure"
    )
    assert "callmpi_allreduce(" in guard and "mpi_max" in guard, (
        f"{label} allocation failure must reach collective consensus before payload use"
    )
    assert "global_bad" in guard and "return" in guard, (
        f"{label} collective allocation failure must return before payload use"
    )


complement_body = subroutine_body(
    COMPLEMENT.read_text(errors="replace"),
    "materialize_dg_hybrid_projected_pw_tile",
)
require_checked_local_allocation(complement_body)

stream_body = subroutine_body(
    STREAM.read_text(errors="replace"),
    "initialize_dg_hybrid_fragment_basis_stream",
)
require_collective_allocation_phase(
    stream_body,
    "fragment_presence(fragment_count)",
    "fragment_presence=0",
    "fragment-presence workspace",
)
require_collective_allocation_phase(
    stream_body,
    "basis%global_ids(nlocal)",
    "basis%fragment_id=fragment_id",
    "stream payload",
)
require_collective_allocation_phase(
    stream_body,
    "ownership(size(wannier_owner)+size(pw_owner))",
    "ownership=0",
    "basis-ownership workspace",
)

print("Task 5 allocation safety source contract: PASS")
