# Portable GS Archive v2 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement an interoperable Zarr-v3 GS array contract, corruption checks and Fortran stream interchange without collecting actual GS data.

**Architecture:** A new explicit `salmonpy.gs_archive` component writes/opens schema `salmonpy.gs/2`. Periodic unpolarized real-space arrays use little-endian float64, wavefunctions split into real/imaginary arrays with axes(k,band,z,y,x), chunks(1,orbital_block,nz,ny,nx). A final manifest indexes exact layouts and SHA256 of every file. Scalar physical metadata, structures, kpoints, weights, occupations, density and potential references have explicit units/conventions. Lazy block reads recheck hashes; no full-wavefunction allocation on open. User budgets are optional bytes for files/decoded blocks, not machine-wide memory limits. Existing native GS/1 and TDDFT paths remain until an audited runtime converter is integrated.

**Tech Stack:** NumPy, stdlib JSON/hashlib/pathlib, strict Zarr-v3 regular chunks/bytes codec, ISO_FORTRAN_ENV stream helper.

## Approved scope

User approved the format-first approach, synthetic round-trip/corruption tests then investigation of SALMON-side connection. Do not disguise native compiler-specific records as portable arrays. Inspect native I/O and implement a standalone Fortran chunk adapter; report native GS-result export/import as unsupported until actual SALMON orchestration and physical validation are available. No full SALMON source import, GS collection, remote job changes or binary replacement. The new archive is array interchange, not a complete SCF/RT-history checkpoint or automatic TDDFT restart.

## Tasks

1. Fail-first tests for bitwise round-trip, edge padding, random block reads, physical metadata, potential references only, corruption/traversal/symlink/extra-file rejection, finite-value checks, optional budgets and failed-write cleanup.
2. Implement strict metadata and array validation, one-k/one-block streaming writer, manifest-last publication, lazy verified reader. Preserve input arrays; readonly returned metadata/arrays. Register public API and package guides. Validate against independent Zarr reader when available.
3. Inspect pinned native read/write_wavefunction and info/occupation formats. Add standalone standard Fortran little-endian real64 chunk read/write helper with scalar IEEE validation and checked allocations. Compile a synthetic interoperability probe only after source audit; audit every app object before linking. Test both directions and rejected NaN/Inf/wrong-size/overwrite.
4. Full source and installed-wheel suites; source-distribution→wheel retains native prohibition checks. Save only synthetic metadata/results/probe logs and distribution artifacts in canonical calc2026. Review and commit; document remaining native GS conversion/TDDFT linkage explicitly.

## Completed and remaining

The format-first component, lazy k/band reads, immutable metadata, explicit budgets, references-only potentials, public catalog/guides and standalone Fortran raw-byte helper are implemented.255 source/installed tests and independent Zarr3.4.2 parity validated; native Fortran probe passes native/swap/big_endian runtime I/O settings, NaN/Inf/size/overwrite rejection. Reviewer fixed reduced-coordinate axis semantics, manifest-first early unknown-tree rejection and runtime-endian conversion sensitivity; final review found no blockers. All data are synthetic, no GS collected.

Native SALMON checkpoint code was inspected: wfn.bin has real/complex and MPI-specific stream layouts, info/occupation use sequential unformatted records. Full GS/1↔GS/2 runtime adapters, physical-material validation and GS/2 TDDFT initialization remain unimplemented. Do not call the current component a complete native restart migration. Canonical records: calc2026/salmonpy/2026-11-14-gs-archive-v2.
