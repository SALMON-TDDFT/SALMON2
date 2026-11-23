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

## Completed milestone (2026-11-15)

- Fast PR c9e8ce43 source audit passed; the PR checkout remained clean/unchanged. New external reference/overlaid/final snapshots preserve mandatory configure/incremental source and post-link application-object audits. No driver prohibition bypass or upstream merge.
- New-type native probe failed first for missing module. Corrupt stored PT-CN dt initially produced STOP34; scalar validation was added and preserves state on rejection.1/2OMP probes pass.
-35private ACE/history/cache fields moved under an explicit type; two lifecycle operations accept it as an argument. Legacy adapter still uses one default instance. Geometry/kernel/Wannier/public flags and global configuration remain future extraction work.
- Actual pinned ACE/history algorithms validate independent factors, interleaved accepted histories/prediction, two dt values, both destroy/recreate orders, rejected stage/dt atomicity, full-source lifecycle and exact allocated-payload accounting. All allocated variants total4616bytes on GNU, including history-array descriptor/fixed elements; clear returns0. This is not RSS.
- Whole HSE solver:4regular k meshes ×2dt ×4modes (Taylor/off, Taylor/learned, direct full, learned PT-CN),12steps each.1/2OMP base vs changed current/energy formatted output files and scientific factor/residual traces match. This does not certify internal-bit or physical GS accuracy. Learned Taylor predictions are actually installed (24total steps across the8learned cases); regression now requires nonzero prediction use and PT-CN residual traces.
- Source292tests and final sdist-built installed-wheel292tests passed, no skips. Installed prep/probe and32case2OMP whole-solver repeat pass. Wall-time records are end-to-end diagnostics, not a statistically controlled speed claim.
- Failed smoke fixtures are preserved: explicit file_kw invalidates the PR full regular-mesh assumption; builtin regular meshes fix input without relaxing validation. Learned PT-CN also requires explicit serial parallel-decomposition input before its admission gate, not automatic zero defaults. No real GS calculation/collection was needed.
- Independent review found no numerical transformation/ownership blocker. Exact accounting and intended-fast-path execution checks were added from feedback. New overlay/experimental helper remains separate from the default SALMONpy runtime, and FDTD continues common dt.
- No production jobs/fixed binaries/PR checkout were changed. No full solver source import, push or performance/whole-native-multi-material claim.
