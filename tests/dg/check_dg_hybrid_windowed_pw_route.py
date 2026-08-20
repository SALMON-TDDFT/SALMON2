#!/usr/bin/env python3
"""Source contract for the narrow windowed-PW/Wannier hybrid route."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TYPE_SOURCE = ROOT / "src/common/dg_hybrid_windowed_pw_types.f90"
BASIS_SOURCE = ROOT / "src/common/dg_hybrid_windowed_pw_basis.f90"
SELECTION_SOURCE = ROOT / "src/gs/dc/dg_hybrid_wannier_selection.f90"
COMPLEMENT_SOURCE = ROOT / "src/common/dg_hybrid_wannier_complement.f90"
COMMON_CMAKE = ROOT / "src/common/CMakeLists.txt"
DC_CMAKE = ROOT / "src/gs/dc/CMakeLists.txt"

assert TYPE_SOURCE.is_file(), "missing windowed-PW hybrid catalog types"
assert BASIS_SOURCE.is_file(), "missing covariant windowed-PW basis primitive"
assert SELECTION_SOURCE.is_file(), "missing Wannier symmetry-block selector"
assert COMPLEMENT_SOURCE.is_file(), "missing local Wannier orthogonal-complement primitive"

type_source = TYPE_SOURCE.read_text().lower()
cmake_source = COMMON_CMAKE.read_text().lower()
dc_cmake_source = DC_CMAKE.read_text().lower()

for required_contract in (
    "module dg_hybrid_windowed_pw_types",
    "type,public::s_dg_hybrid_pw_packet",
    "type,public::s_dg_hybrid_basis_catalog",
    "accepted_wannier_blocks",
    "rejected_wannier_blocks",
    "wannier_fingerprint",
    "window_fingerprint",
    "packet_fingerprint",
    "catalog_fingerprint",
):
    assert required_contract.replace(" ", "") in type_source.replace(" ", ""), (
        f"missing hybrid catalog contract: {required_contract}"
    )

assert "dg_hybrid_windowed_pw_types.f90" in cmake_source, (
    "windowed-PW hybrid catalog must be part of the SALMON build"
)
assert "dg_hybrid_windowed_pw_basis.f90" in cmake_source
assert "dg_hybrid_wannier_complement.f90" in cmake_source
assert "dg_hybrid_wannier_selection.f90" in dc_cmake_source

production_source = "\n".join(
    path.read_text(errors="replace").lower()
    for path in (ROOT / "src").rglob("*.f90")
)
for obsolete_global_contract in (
    "s_dg_wpw_s_orthogonal_complement",
    "a_owned_w_global_p",
    "dg_wpw_production_context",
):
    assert obsolete_global_contract not in production_source, (
        "the new hybrid route must not restore obsolete global WPW state: "
        f"{obsolete_global_contract}"
    )

print("windowed-PW hybrid route contract: PASS")
