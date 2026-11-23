# TDDFT replay snapshots implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Expose in-memory snapshot/restore for process-isolated TDDFT and verify fixed-dt continuation.

**Architecture:** Reuse immutable StateSnapshot tokens and the common compatibility/integrity validator. Restore starts an independent worker from the state's owned initial bundle, replays accepted A history, checks every endpoint current, and commits by exchanging worker resources only after validation. Retirement errors must not turn an already committed restore into a reported restore failure.

**Tech Stack:** Python, NumPy, subprocess worker protocol, unittest; separately audited unchanged SALMON runtime.

## Task 1: Public contract tests

Create `salmonpy/tests/test_tddft_snapshots.py`. Use a real isolated deterministic Python protocol worker (portable launch via Python on Windows). Test initial rewind, divergent branches, 20-step continuation, token immutability/validation, same model and dt, closed state, spawn/replay/cleanup failure atomicity, and finish output history. Run `python3 -m unittest salmonpy.tests.test_tddft_snapshots`; verify missing API failures before implementation.

## Task 2: Replay transaction

Modify `salmonpy/tddft.py`. Retain initial model identity, snapshot A/current histories, fingerprint fixed settings, validate variable-length records through `_SnapshotMixin`, replay in a separate worker, then swap resource ownership. Keep replayed currents. Retire workers with fallback cleanup and retain private cleanup diagnostics without raising after commit. Preserve existing failed-step recovery and finish behavior.

## Task 3: Native continuation verification

Create `tools/test_tddft_snapshots.py`. Load an unchanged runtime through checksum verification, use an existing synthetic-orbital input bundle, compare 20 subsequent currents bitwise after same-state and independent-state restoration, rewind to initial state, and reject a different dt. Save only current traces and metadata to an external calculation-results directory. Do not collect GS data or build a solver.

## Task 4: Documentation and regression

Update Japanese/English SNAPSHOTS guides, AI_USAGE, and matching manual.json embedded guides. Explain replay cost, two simultaneous workers, fixed settings, current tolerance, branch output regeneration, and exclusions (persistent disk/embedded/MPI restore). Run `python3 -m unittest discover -s salmonpy/tests` with checksum-verified existing kernels; review diff and run `git diff --check`.
