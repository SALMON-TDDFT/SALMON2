#!/usr/bin/env python3
"""Protect current bounded full-cell projection lifetimes.

The superseded test expected move_alloc(global_seed_values,w90_anchors), a
symbol absent from the bridge snapshot.  The reachable implementation bounds
the full-cell Hpsi work by tiles and releases transient component storage.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()


def subroutine(source: str, name: str) -> str:
    start = source.index(f"subroutine {name}")
    end = source.index("end subroutine", start)
    return source[start:end]


def violations(source: str) -> list[str]:
    body = subroutine(source, "ow_build_hamiltonian")
    compact = "".join(body.replace("&", "").split())
    problems: list[str] = []
    if body.count("call project_dg_full_cell_hamiltonian_tiles") != 3:
        problems.append("full/kinetic/local projections must each use the tiled projector")
    if compact.count("min(16,nwann)") != 3:
        problems.append("every full-cell projection must cap its destination tile width")
    for output in ("full_rows", "kinetic_rows", "local_rows"):
        if output not in body:
            problems.append(f"missing local projected component {output}")
    sym = body.find("call symmetrize_dg_distributed_pencil_rows")
    release = body.find("deallocate(component_rows)")
    if sym < 0 or release < sym:
        problems.append("component_rows must be released immediately after symmetrization")
    for name in ("kinetic", "local", "nonlocal"):
        old_release = body.find(f"if(allocated(ow_last_{name}_rows))deallocate(ow_last_{name}_rows)")
        replacement = body.find(f"allocate(ow_last_{name}_rows,source={name}_rows)")
        if old_release < 0 or replacement < old_release:
            problems.append(f"ow_last_{name}_rows replacement does not release its previous owner")
    if "allocate(global_h" in body or "global_h(:,:)" in body:
        problems.append("Hamiltonian callback must not own a replicated global dense matrix")
    return problems


assert not violations(SOURCE), violations(SOURCE)

# Mutation fixtures prove that the contract observes both the memory bound and
# the transient-lifetime release rather than merely finding the routine name.
unbounded = SOURCE.replace("min(16,nwann)", "nwann", 1)
assert any("cap its destination tile width" in item for item in violations(unbounded))
leaked = SOURCE.replace("deallocate(component_rows)", "! removed by mutation fixture", 1)
assert any("component_rows must be released" in item for item in violations(leaked))

print("PASS current bounded full-cell projection memory lifetimes and mutation fixtures")
