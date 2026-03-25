#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


CHECKS = [
    {
        "name": "task1_inventory",
        "patterns": [
            {
                "file": "src/rt/dg/rt_dg_integrator_derivative.f90",
                "regex": r"dg_frag%(?:S_mat(?:_prop)?|H_mat|momentum_mat)(?:_c)?\(",
                "message": "dense fragment-only operator access remains in derivative path",
            },
            {
                "file": "src/rt/dg/rt_dg_integrator_unitarity.f90",
                "regex": r"dg_frag%(?:S_mat|S_mat_prop|S_mat_c|S_mat_prop_c)\(",
                "message": "dense fragment-only overlap access remains in unitarity path",
            },
            {
                "file": "src/rt/dg/rt_dg_basis_projection.f90",
                "regex": r"dg_frag%(?:H_mat|S_mat(?:_prop)?|momentum_mat)(?:_c)?\(",
                "message": "dense fragment-only operator access remains in basis projection",
            },
            {
                "file": "src/rt/dg/rt_dg_fragment_basis_update.f90",
                "regex": r"dg_frag%(?:H_mat|S_mat(?:_prop)?|momentum_mat)(?:_c)?\(",
                "message": "dense fragment-only operator access remains in basis update",
            },
            {
                "file": "src/rt/dg/rt_dg_fragment_basis_update_soi.f90",
                "regex": r"dg_frag%(?:H_mat|S_mat(?:_prop)?|momentum_mat)(?:_c)?\(",
                "message": "dense fragment-only operator access remains in SOI basis update",
            },
        ],
    },
    {
        "name": "task2_overlap_broadcast",
        "patterns": [
            {
                "file": "src/rt/dg/rt_dg_fragment_ops.f90",
                "regex": r"subroutine\s+ensure_overlap_prop_available[\s\S]*?(?:comm_bcast|sync_dense_matrix_to_blocks_runtime)[\s\S]*?end\s+subroutine\s+ensure_overlap_prop_available",
                "message": "steady-state overlap helper still broadcasts or rebuilds persistent dense overlap state",
            },
            {
                "file": "src/rt/dg/rt_dg_fragment_hamiltonian.f90",
                "regex": r"overlap_prop_root_authoritative\s*=\s*\.true\.",
                "message": "non-SOI Hamiltonian still marks dense overlap state as root-authoritative",
            },
            {
                "file": "src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90",
                "regex": r"overlap_prop_root_authoritative\s*=\s*\.true\.",
                "message": "SOI Hamiltonian still marks dense overlap state as root-authoritative",
            },
        ],
    },
    {
        "name": "task3_integrators",
        "patterns": [
            {
                "file": "src/rt/dg/rt_dg_integrator_unitarity.f90",
                "regex": r"matmul\s*\(\s*(?:cmplx\s*\()?\s*dg_frag%(?:S_mat|S_mat_prop|S_mat_c|S_mat_prop_c)",
                "message": "unitarity still applies dense fragment-only overlap directly",
            },
            {
                "file": "src/rt/dg/rt_dg_integrator_derivative.f90",
                "regex": r"S_eval\(:,\s*:\)\s*=\s*(?:cmplx\()?dg_frag%(?:S_mat|S_mat_prop|S_mat_c|S_mat_prop_c)",
                "message": "derivative path still copies dense fragment-only overlap directly",
            },
        ],
    },
]


def iter_matches(text: str, regex: str) -> list[tuple[int, str]]:
    pattern = re.compile(regex, re.MULTILINE)
    matches = []
    for match in pattern.finditer(text):
        line_no = text.count("\n", 0, match.start()) + 1
        line = text.splitlines()[line_no - 1].strip()
        matches.append((line_no, line))
    return matches


def main(argv: list[str]) -> int:
    enabled = set(argv[1:]) if len(argv) > 1 else {check["name"] for check in CHECKS}
    findings = []

    for check in CHECKS:
        if check["name"] not in enabled:
            continue
        for item in check["patterns"]:
            path = ROOT / item["file"]
            text = path.read_text()
            matches = iter_matches(text, item["regex"])
            for line_no, line in matches:
                findings.append(f"{check['name']}: {item['file']}:{line_no}: {item['message']}: {line}")

    if findings:
        sys.stderr.write("DG fragment dense elimination audit failed:\n")
        sys.stderr.write("\n".join(findings) + "\n")
        return 1

    print("DG fragment dense elimination audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
