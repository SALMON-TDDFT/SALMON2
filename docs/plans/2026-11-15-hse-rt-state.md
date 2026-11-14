# Fast-PR HSE RT state extraction implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extract ACE/factor-prediction/Taylor/PT-CN mutable RT state from the user's fast PR revision into an explicitly owned type and lifecycle operations, without changing exchange numerics.

**Architecture:** Pin local fast PR c9e8ce431f43eefd4a88998134adb4796f8ec837. An opt-in external overlay adds exx_rt_state, moves selected private ACE/history/workspace globals into s_exx_rt_state and makes the legacy adapter use one default instance. Lifecycle routines take an explicit state/config argument and return status before mutation on invalid stage/dt. Existing exchange loop/FFT/communication code and public legacy entry points remain. This is a transitional extraction: default native adapter, geometry/kernel/Wannier state, public legacy flags and configuration remain shared; no full multi-material claim, global swap, opaque-handle copy or native ABI1 claim.

**Tech Stack:** external pinned PR archive, small opt-in patch and module, source/object prohibition checks retained in fresh CMake builds, actual exx_ace/exx_factor_history native tests and standalone/full-solver regression where available. No full source import, PR edits, jobs, physical GS collection or push.

## Steps

1. Write a native fail-first test requiring the new type/helpers, using actual pinned ACE/history modules. Interleave two states with different factors/history/dt; compare each with its own independent reference, Taylor trial acceptance/reset, PT-CN constant-dt rejection and sibling destruction. Run old snapshot and verify missing module/API failure.
2. Implement state container/lifecycle functions. Replace only selected module-state references in exx_native; preserve public wrappers and numerical order. Keep persistent shared resources explicitly outside this first extraction. Add upstream attribution and pinned source hashes.
3. Prepare a fresh external PR snapshot, audit it before any build, retain/reinstate configure and incremental object/source guards locally. Never build the PR checkout in place or disable mandatory audits. Audit every application object before ad-hoc native probe links.
4. Compare native factor/history/ACE results against unchanged routines; validate finite/NaN/Infinity, memory/lifetime behavior and failure atomicity for helper validation. Build full HSE solver and run available synthetic non-GS regression in new directories; distinguish compile-only, lifecycle proofs and actual exchange accuracy/performance. Do not invent a speed claim.
5. Source/installed packaging tests, review, canonical calc2026 records and local commit. Document remaining geometry/config/public flag/worker boundaries before full native multi-context support.
