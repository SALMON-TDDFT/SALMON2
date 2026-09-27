# Bounded DC-to-mesh reconstruction implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove global-grid orbital and coverage scratch from complex DC-to-mesh initialization.

**Architecture:** Fragment owners retain existing validated fragment input. Broadcast each destination's grid/orbital bounds (eight integers). Validate core coverage and contract fragment coefficients in flat local-grid chunks of at most65536 points. Reduce each chunk only to its destination in icomm_ro. Never allocate a global coverage or orbital array. Preserve serialized metadata checks, fragment ownership and projected LCFO configuration. This changes temporary reconstruction memory; resident fragment basis and coefficient data are not streamed yet.

**Tech Stack:** Fortran pure tile kernels, existing rooted comm_summation, Python/Fortran and native MPI regressions.

1. Add failing tile-kernel probe: complex coefficients, wrapped/nonmonotone maps, disjoint tiles, chunk crossing row/plane boundaries, missing/duplicate coverage and zero-band fragments. Add native scratch diagnostic expectation to existing MPI1/2/4 pulse parity test and observe failure before code.
2. Add src/gs/dc/lcfo_mesh_tile.f90 and CMake entry. Extract bounded coverage and contraction. In lcfo_complex.f90 remove full-grid allocations; validate all files first, then chunk coverage, load fragment records as before, configure projected basis unchanged, reduce reconstructed chunks to grid/orbital owners.
3. Verify tile oracle against independent dense reconstruction; run native spatial pulse/work/water, existing projected LCFO and malformed-file/provenance regressions; HSE ON/OFF builds and independent review. Save scratch sizes and parity, document remaining fragment-resident memory and local commit.


## Execution record

- Observed unit failure for missing tile module and native failure for absent bounded-scratch diagnostic before implementation.
- Implemented bounded pure tile kernels and destination-rooted reconstruction/coverage reductions; metadata and projected configuration preserved.
- Tile oracle, 15 Ehrenfest tests, 8 input tests, HSE ON/OFF builds and existing projected LCFO response/provenance/orbital-layout suite passed.
- Saved pre-change trajectory comparison: H4 MPI1/2/4 unchanged at output precision; water MPI2 energy change1.03e-12 eV.
- Independent review found no blockers. Ruling: defer repeated chunk intersection scans; this stage bounds memory and does not claim improved startup time. Resident fragment input is unchanged.
- Evidence: docs/results/pbeh40-rvv10/dc-tile-validation.json. No push or merge.
