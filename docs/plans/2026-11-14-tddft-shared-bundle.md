# TDDFT Shared Replay Bundle Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reduce managed initial GS bundles from 2N+1 to N+1 for N process-isolated states, preserving replay, independent workers and checksums.

**Architecture:** The initial loader owns one private baseline store. Each state leases this store and owns a separate writable worker copy. Closing the loader refuses new states but retained leases preserve replay; the baseline is deleted after the loader and final state release it. Restore candidates fork an existing state lease, including after loader close. Every spawn copies to an independent working directory and compares its full manifest before launching a process. No symlinks/hardlinks or filesystem-specific clone behavior.

**Tech Stack:** Python tempfile/threading ownership, existing process/replay/checkpoint APIs and SHA256 manifest.

## Approved design and scope

User approved the shared immutable baseline / independent worker design on continuing after design review. Reflinks are deferred for portability; keeping independent state baselines would retain duplicate disk I/O. No SALMON native changes/builds, physical calculation, GS collection or remote job changes. Synthetic protocol fixtures and generated opaque payloads test copies/lifecycle, not SALMON physics.

## Tasks

1. Add failing tests: shared baseline paths but independent worker inodes; close initial/refuse new states; continued stepping/replay; final lease cleanup; spawn/candidate failure reference cleanup; corrupted baseline rejected before process spawn; changed worker files do not alter baseline/other workers.
2. Add private store/lease classes to tddft.py, thread-safe acquisition/release and idempotent cleanup; acquire in state creation, fork in restore; use existing stop/retirement paths. Verify copied manifests in _spawn and remove extra baseline copies. Preserve atomic restoration and retained-output recovery lifetimes.
3. Run full source suite against unchanged hash-validated current LD/SBE library, protocol replay/checkpoint/output tests, and review. No new physical SALMON build/execution.
4. Benchmark literal prechange vs new implementation with deterministic workers and generated non-GS payloads, 1/4/8 states, independent subprocess repeats; logical copy bytes/physical allocation/time/RSS reported separately, numerical parity required. Update bilingual guide and manual embeddings. Rebuild sdist→wheel with mandatory native audits, installed suite outside checkout; save canonical records, commit without push.

## Completed validation

246 source and installed-wheel tests pass. Shared baseline lifetime, failure/refcount cleanup, GC, worker-file isolation, prelaunch corruption rejection, replay after loader close and retained-output recovery covered. With16MiB opaque payloads and1/4/8 deterministic protocol workers, all9 old/new comparisons are bitwise identical. For8states, managed payload copies17→9 (272→144MiB), state creation0.2570→0.1825s and restore0.3104→0.2322s; these synthetic results do not claim physical TDDFT speedups or physical filesystem extent savings. No GS data collected, native/runtime sources unchanged and existing binaries/jobs untouched. Independent review found no blockers. Canonical records: calc2026/salmonpy/2026-11-14-tddft-shared-bundle.
