# Native Gamma GS to GS/2 Export Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Export a verified serial Γ GroundStateResult to portable GS/2, preserving wavefunctions and occupations, recording reconstructed density provenance.

**Architecture:** Extend the explicit compiler-paired converter build with a native record exporter supporting real/complex one-rank Γ streams. The exporter reads native records with the same GNU profile and outputs portable little-endian streams. Python validates source manifest/physical input/eigenvalue AU format, copies verified records into private temporary input, maps output streams, reconstructs density in fixed band order, and writes GS/2 by orbital blocks. Optional typed export provenance records source hashes, native orbital kind, local channels and derived density semantics. GroundStateResult.export_archive is the public entry; no actual GS data collected during tests.

**Tech Stack:** GNU Fortran checked scalar IEEE and raw bytes, existing audited build tool/native writer probes, NumPy memmap and GS/2 writer.

## Approved scope and limits

User approved Γ/serial real+complex records, AU eigenvalue parsing, density reconstructed from u/occupation, and source metadata checks. Density mixing history is not migrated; this is not exact SCF-state restart. Only DFT-generated GS/1 inputs supported; reject inconsistent/missing metadata/records, other k/spin/MPI/units. No guess/fill of missing eigenvalues. Existing binaries/jobs and shared kernel/runtime sources unchanged. Synthetic complete-result fixtures do not establish real-material convergence.

## Tasks

1. Fail-first tests for public method, source/destination/metadata/units/provenance and cleanup; preserve backward compatibility for old converter import-only records and existing GS/2 files.
2. Add Fortran exporter with strict header/record sizes/kind/finiteness/capacity checks and little-endian output. Extend explicit builder to audit sources/all objects before links and test unchanged native real/complex writer→exporter→converter/read parity, invalid files and runtime byte-swap settings.
3. Implement Python export with source input/model checks, exact pseudo/eigen metadata, bounded memmap/band reconstruction, immutable provenance; no full-wavefunction read. Handle Windows mapping cleanup before temp deletion. Add public GroundStateResult entry and documentation.
4. Synthetic end-to-end native result→GS/2→native TDDFT, currents/restore and real+complex byte parity; source and installed-wheel suites, installed builder. Save canonical synthetic-only records and commit without push.

## Completed validation

272 source/installed tests pass. Explicit builder verifies unchanged native real/complex writer→export portable values and occupations; old import-only converter records remain accepted for import, export requires new probed exporter. Synthetic GS/1 real+complex→GS/2→native actual unchanged-runtime six-step currents and restored branch match bitwise; native/swap I/O variants and source/installed variants agree. Derived density integral4 electrons within binary64 rounding; energy/eigenvalues not recalculated. Source GS files and fixed native/runtime binaries unchanged. No actual SCF/GS collection or real-material convergence claim.

Independent review found no remaining blockers. Source archive preservation, metadata/input/hash/unit/occupation mismatch, temporary cleanup and typed provenance covered. Canonical records: calc2026/salmonpy/2026-11-15-gs-gamma-export. Multi-k/spin/MPI/foreign compiler profiles and physical SCF validation remain out of scope.
