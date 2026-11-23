# Serial TDDFT shared-configuration session implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Validate multiple electronic states inside one native process, with one immutable material/grid/XC/dt configuration and independently owned mutable state.

**Architecture:** A separate experimental serial session library initializes a baseline once and creates native state handles from that baseline. The first implementation uses builtin PZ to avoid copying opaque Libxc handles. Context arrays, iteration parity, A/current and scratch output units are private; shared PZ workspace and configuration remain session-owned and calls are serialized on the main thread. One native session lifetime per process remains explicit; no MPI, concurrent threads, differing materials/dt, native checkpoints or rollback-on-native-fatal-error claim. Existing process API is unchanged and remains the reference.

**Tech Stack:** Existing external fixed SALMON source + compatibility patch, a separately distributed small adapter, explicit audited fresh GNU/PIC build, ctypes experimental Python wrapper. No full solver import or production changes.

## Superseded scope

The user subsequently selected design for different materials and dt in one process. See `2026-11-15-tddft-independent-contexts-design.md`. The local same-configuration prototype was unvalidated, its build stopped, and its files retained outside the checkout; no public feature was delivered. The tasks below are historical planning, not completion claims.

## Original accepted design

User approved the same material/grid/XC/dt scope after review of the context audit. Work proceeds inline in the existing worktree. Start with PZ and prove ownership before extending XC. Keep the baseline (and thus N+1 electronic copies) until session close; do not optimize sharing mutable substructures prematurely. Initialize legacy input only once; never shallow-copy external XC resources. Separate scratch units avoid legacy output collisions. Session close refuses live states; closing/recreating states does not reinitialize global configuration. Handles are monotonic 64-bit IDs and stale handles never refer to a replacement.

## Steps

1. Add failing wrapper contract tests: capability/hash/scope, thread/parameter rejection, lifetime, sibling independence and stale handles. Confirm missing module/API failure.
2. Add native session adapter and explicit builder. Fresh source snapshot only; source audit before configure/build and application object audit before ad-hoc link; retain CMake audits. Runtime/shared kernels remain unchanged. Native step validation is scalar IEEE and happens before state mutation.
3. Implement experimental Python Session/state wrapper using an immutable owned copy of TDDFTInitialState; shared dt/capacity; finite A, thread and library checks; no import-time native load/build. Restrict builtin PZ inputs. Failed validation leaves handles/current/time unchanged. Native fatal errors remain experimental process risk.
4. Actual native tests for2/3/6/8 and weighted2 synthetic GS fixtures, two different interleaved histories versus unchanged independent worker processes; close either state, recreate, exhaust capacity, reject stale IDs/NaN/stale A. Measure RSS/lifetimes without universal hard cap. Run each mesh in a separate host process because session global configuration has one lifetime.
5. Run source/installed tests, sdist→wheel with audits, package tool/adapters/docs. Save synthetic-only records in canonical calc2026, review numerical/lifetime limits, commit locally without push. Existing jobs/fixed executable untouched.
