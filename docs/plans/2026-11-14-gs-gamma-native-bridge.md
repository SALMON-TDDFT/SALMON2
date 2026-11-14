# GS/2 Gamma Native Bridge Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Connect portable GS/2 Γ-point arrays to a one-rank unpolarized native SALMON TDDFT initialization without changing fixed runtimes.

**Architecture:** An explicit audited converter build derives GNU compiler and active flags from the target runtime's CMake cache, binds its executable SHA256, and writes compiler-native info/occupation records plus a complex wavefunction stream. A standalone probe compiles unchanged native read_wavefunction/write_wavefunction procedures from an audited external source tree with minimal type scaffolding and compares converter output bytewise with the native writer. Python verifies target/converter hashes and potential hashes, streams GS/2 bands into converter inputs, creates a new independent RT bundle, and returns TDDFTInitialState. Conversion I/O environment is recorded and checked before each later spawn. No compiler-format guessing in Python.

**Tech Stack:** Python/NumPy, trusted runtime build record/CMake cache, standard Fortran stream adapter, existing immutable bundles and process TDDFT API.

## Approved scope

User approved Γ point, one rank, periodic unpolarized first; other k/spin/MPI layouts are rejected. Runtime and converter must share GNU compiler/profile in this initial implementation. No implicit build/download, physical GS collection, fixed executable replacement, remote-job mutation, or full SALMON source import. Synthetic native read/write/control validation is distinguished from SCF/physical convergence. GS/1→GS/2 export remains subsequent work.

## Tasks

1. Fail-first public entry and input/hash/unsupported-scope tests. Preserve baseline and existing binary hashes.
2. Add explicit converter build tool, Fortran converter and synthetic probe scaffolding. Audit external source tree before extracting native readers/writers; audit snapshotted sources before compiling, all application objects before linking. Refuse historical sources with prohibited references; preserve evidence.
3. Implement Python GammaRestartConverter loader, deterministic RT bundle preparation and GroundStateArchive.tddft_initial_state. Validate Γ/exactweight, grid/electron/XC scope, pseudopotentials and target SHA before conversion. Preserve bounded band I/O, no source archive mutation, no overwrite and failed-bundle cleanup. Bind environment receipt and enforce at every spawn.
4. Test unchanged native reader/writer byte parity and synthetic full-runtime initialization/steps when feasible, plus snapshot restore. Run source/wheel suites, installed-tool build/packaging checks. Document supported scope and limitations, save canonical synthetic records and commit without push.

## Completed validation

265 source and installed-wheel tests pass. Source and installed converter builders derive the paired GNU runtime profile, audit safe full source/snapshots/all objects, and pass unchanged native wavefunction reader/writer byte parity. Real unchanged SALMON runtime initialized synthetic orthonormal Γ orbitals; six-step independent/restored branches and source/installed/native/swapped-I/O variants have identical currents. These are synthetic connection tests, not SCF or physical GS accuracy validation. Existing fixed runtime/native kernel binaries and shared runtime patch unchanged.

Review fixed16-byte complex output overflow checks, capturedenv races, pre-mkdir environment validation, and full-wfn hash allocation; added regressions. Portable source archive is retained as checkpoint-protected input. Pseudopotential local-channel/XC/original-GS consistency remains caller responsibility. Historical embedded-base tree failed source prohibition audit and was neither built nor executed; safe force-source passed. Failed diagnostic link attempt/source/logs preserved. Canonical results: calc2026/salmonpy/2026-11-14-gs-gamma-native-bridge.

Native GS/1→GS/2 export, multiple k points/spins, MPI/native compiler profiles other than tested GNU, and converged-material validation remain separate work.
