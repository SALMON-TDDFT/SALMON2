# TDDFT disk checkpoints implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Persist accepted process-isolated TDDFT replay history and restore it in a fresh process using an explicitly supplied identical initial model.

**Architecture:** `TDDFTState.save_checkpoint(directory)` writes a new Zarr v3 group with uncompressed little-endian float64 arrays and a final hash manifest. `restore_checkpoint(directory, max_bytes=None)` verifies all metadata, capacity-derived byte ceilings and optional caller budgets, exact file inventory, checksums, model/input/initial-file hashes, dt, capacity, threads and executable hash before handing a rebound immutable token to existing transactional replay. No executable paths are loaded or executed from the checkpoint. Initial GS and pseudopotential files are not copied; callers supply their existing initial bundle. In-memory tokens retain their stricter same-object contract.

**Tech Stack:** Python standard library, NumPy, existing TDDFT worker/replay; no required Zarr dependency.

## Accepted scope and tradeoffs

The user approved disk persistence after explicitly deprioritizing structural optimization. Native wavefunction checkpoints would avoid replay but require substantially more native state work. Saving history plus a reference fingerprint reuses verified replay and is selected for this first implementation. Persistence of light models and a self-contained GS archive remain separate features. Input/model settings are content-bound rather than process identity-bound; runtime SHA and threads must match exactly. Same-directory overwrite is rejected. Interrupted writes lack a valid final manifest; ordinary failures clean their new destination. Power-loss durability and concurrent writers are not promised.

Zarr metadata follows https://zarr-specs.readthedocs.io/en/latest/v3/core/ and the bytes codec specification. History is chunked along time in blocks of at most1024 vectors, padded to the full chunk extent. An empty history has no chunks. Metadata and all chunks are inventoried and hashed; these hashes detect corruption, not malicious publisher identity.

## Task1: Public contract and failure tests

Create `salmonpy/tests/test_tddft_checkpoints.py`, reusing the deterministic protocol-worker fixture. Test fresh model/process continuation, zero-step checkpoints, chunk boundaries, corruption, dt/model/executable/thread incompatibility, bounded loads, overwrite refusal, interrupted save cleanup, replay failure atomicity, path/symlink rejection, and absence of GS payloads. Run before implementation and confirm missing-API failures.

## Task2: Minimal persistence adapter

Create `salmonpy/_tddft_checkpoint.py`; add two forwarding methods in `salmonpy/tddft.py`. Keep archive validation before any worker spawn. Construct a StateSnapshot with the explicitly supplied model owner only after all integrity checks; reuse existing restore for replay and atomic commit. Do not modify native sources, runtime patch or existing jobs.

## Task3: Documentation and distribution

Update paired SNAPSHOTS guides, API catalog/manual and AI guide; install the private adapter in wheel CMake. Extend native snapshot tool to save, close/recreate model, restore from disk and compare20 future steps, retaining only history and scalar results. Run targeted and full Python suites. Build sdist→wheel with retained prohibition audits and run installed-package tests outside checkout. Verify compatibility with an independent Zarr reader where available.

## Task4: Verification and records

Use an unchanged previously audited runtime and synthetic input bundle; verify recorded SHA and compare continuation currents bitwise. No GS data collection. Save results under canonical calc2026, distinguish interface regression from physical validation. Review diffs, commit completed changes, do not push or modify existing binaries/jobs.

## User clarification on limits

The user challenged arbitrary64MiB and inventory/depth caps. Removed them. Defaults now follow the explicit receiving model capacity,48bytes per accepted step, chunk padding and all metadata/hash overhead. A caller may provide max_bytes as an additional limit. Canonical file/directory membership replaces unrelated entry/depth caps, and manifest file counts must match the format before any expected inventory expansion. These are stored/decoded-payload bounds, not total RSS or automatic permission to consume machine memory.
