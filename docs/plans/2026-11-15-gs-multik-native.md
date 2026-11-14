# Serial Multi-k GS Interchange Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend portable/native GS interchange to aligned orthogonal cells with multiple k points, preserving the existing Gamma API.

**Architecture:** New SerialRestartConverter opt-in requires newly probed multi-k capability; GammaRestartConverter remains exact Gamma. Native converter/exporter add backward-compatible optional nk requests, weighted occupation checks and k-major/band-minor streams. Python writes explicit reduced-coordinate file_kw input, records/verifies its hash, exports strict orthogonal gs_k.data and all eigenvalue blocks, reconstructs weighted density in k/band order. Existing GS/2 axes/chunks unchanged, provenance optionally records k-data source. One serial CPU rank/unpolarized only.

**Tech Stack:** existing GNU compiler-bound audited builder, unchanged extracted SALMON wavefunction readers/writers, NumPy block/mapped I/O and current TDDFT protocol.

## Approved scope

User approved orthogonal serial unpolarized multi-k,2/8-point synthetic tests, wave/occupation bit parity and tolerant current comparison. Orthogonal means native axis-aligned positive diagonal cell (nonorthogonal export uses a different coordinate convention and remains rejected). Existing Gamma public behavior retained. No actual GS collection, SCF/physical certification, native runtime replacement, remote jobs or full source import. GNU paired profile only; no MPI decomposition.

## Tasks

1. Fail-first entry/scope/order/weights/receipt tests. Implement capability-gated subclass and preserve old Gamma converter records/protocol.
2. Add optional nk and weights to native programs with checked16*ng*nb*nk capacity, band×k occupation records, scalar finite/weight checks; extend probe/scaffolding to nk2/8 with distinguishable synthetic orbitals; source/object audits before every new link.
3. Generalize import/export Python paths behind opt-in: exact Γ behavior unchanged; reduced k file, verified receipt, strict orthogonal k export/reciprocal basis/weights, eigen blocks, mapped weighted reconstruction and backward-compatible provenance.
4. Validate unchanged native read/write bit parity, real runtime source→GS/2→native2/8-point currents/restores (rtol1e-12,atol1e-14), nonuniform-weight import and swapped I/O. Full source/installed suites, installed builder and packaging. Review, canonical synthetic-only records, commit without push.

## Completed validation (2026-11-15)

- Source and sdist-built installed wheel: 278 tests each, all passed without skips.
- Audited source and installed converter builds: unchanged native writer/reader parity for Gamma real/complex and multi-k 2/8, including nonuniform occupation/weight proof.
- Actual unchanged serial SALMON runtime, synthetic fixtures: 2/8-point currents bitwise equal; 3/6-point maximum differences 1.734723475976807e-18 / 8.673617379884035e-19, within rtol1e-12 and atol1e-14. Independent nonuniform two-point import bitwise equal. Replay restoration bitwise equal.
- Swapped native record environment and installed-wheel runtime repeat gave identical numerical metrics/currents (build-specific provenance hashes differ).
- No real GS collected; no production jobs or fixed solver binaries changed. General k-coordinate/current bit parity is not promised due text rounding and native weight normalization. GNU/macOS serial scope only; no MPI/Fujitsu/GPU certification.
