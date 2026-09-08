#!/usr/bin/env python3
"""Protect row-owned tiled storage in the current full-cell Hamiltonian path.

The old test required weak-operator row assembly, which the bridge had already
replaced with project_dg_full_cell_hamiltonian_tiles.
"""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()


def subroutine(source: str, name: str) -> str:
    start = source.index(f"subroutine {name}")
    end = source.index("end subroutine", start)
    return source[start:end]


def violations(source: str) -> list[str]:
    builder = subroutine(source, "ow_build_hamiltonian")
    compact = re.sub(r"\s+|&", "", builder)
    problems: list[str] = []
    if builder.count("call project_dg_full_cell_hamiltonian_tiles") != 3:
        problems.append("builder must project full, kinetic, and local row tiles")
    if compact.count("ow_core_values,ow_row_ids,min(16,nwann)") != 3:
        problems.append("each tiled projection must use explicit owned row IDs and a bounded width")
    for declaration in (
        "complex(8),intent(out)::hrows(:,:)",
        "allocate(nonlocal_rows(size(ow_row_ids),nwann)",
        "component_rows(size(ow_row_ids),nwann,3)",
    ):
        if declaration not in compact:
            problems.append(f"missing row-owned extent: {declaration}")
    if re.search(r"(?:allocate\s*\(\s*global_h|global_h\s*\(\s*:\s*,\s*:\s*\))", builder):
        problems.append("builder owns a forbidden replicated global Hamiltonian")
    if re.search(r"allocate\s*\([^\n)]*nwann\s*,\s*nwann", builder):
        problems.append("builder allocates a forbidden Nw-by-Nw local dense matrix")
    if "int(size(hrows),8)*16_8" not in compact:
        problems.append("local-byte evidence must derive from the actual row-owned hrows allocation")
    hermiticity = subroutine(source, "ow_distributed_hermiticity")
    if "row_batch_size" not in hermiticity:
        problems.append("distributed Hermiticity diagnostics must use bounded row batches")
    return problems


assert not violations(SOURCE), violations(SOURCE)

# Mutation fixtures: removing row ownership/bounding or introducing a global
# dense owner must be detected deterministically.
unbounded = SOURCE.replace("ow_core_values,ow_row_ids,min(16,nwann)", "ow_core_values,ow_row_ids,nwann", 1)
assert any("explicit owned row IDs and a bounded width" in item for item in violations(unbounded))
global_dense = SOURCE.replace(
    "logical::finite_t,finite_local,finite_nonlocal,finite_h,update_auxiliary",
    "logical::finite_t,finite_local,finite_nonlocal,finite_h,update_auxiliary\n    complex(8),allocatable::global_h(:,:)",
    1,
)
assert any("replicated global Hamiltonian" in item for item in violations(global_dense))

print("PASS current row-owned tiled Hamiltonian storage and mutation fixtures")
