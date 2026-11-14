# Explicit fast-HSE RT-state routing implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Route the extracted HSE RT state explicitly through exchange/lifecycle call chains; prove that real exx_add_action calls on two supplied states do not use the legacy singleton.

**Architecture:** Add an optional explicit rt_context argument to module entry points and all transitive helpers that access RT state. Each stack frame resolves a local pointer to the supplied state, or to the legacy default when omitted; never swap a process-wide active pointer or copy the whole module. Move remaining public scalar RT flags/energy into the state while retaining legacy aliases to the default instance. Exchange configuration and geometry/kernel/Wannier resources are still shared and remain a subsequent phase; explicit refresh across differing configurations is not yet certified.

**Tech Stack:** New chained opt-in overlay on the pinned c9e8ce43 fast-PR RT-state profile, GNU full-library object-linked native probe, existing32case HSE regression at1/2OMP, source/object audits, wheel validation. No public stable multi-material backend or native C ABI claim.

## Tasks

1. Fail-first native probe calls exx_add_action and lifecycle wrappers with rt_context keyword. Existing module must fail compilation. Use actual compiled HSE objects; audit every application object in both the full build and probe tree before ad-hoc linking.
2. Add explicit routing to all module-level procedures that access RT fields, and propagate it through every internal call. Nested callbacks inherit only their parent's local pointer. Preserve arithmetic/OMP loops and expose no active-global context switch. Preserve legacy default behavior/API for existing callers.
3. Add state-owned public flags/energy with default-only compatibility aliases. Validate two distinct dense ACE operators, interleaved actions, endpoint-averaged trial actions, different PT-CN dt, both close/clear orders, independent buffers and unchanged legacy sentinels. Do not copy opaque FFT/MPI resources; they remain outside this RT-only routing scope.
4. Build a fresh guarded full solver. Compare default32case1/2OMP HSE output/diagnostics to unchanged fast PR. Add installed prep/probe/runtime checks, review, source/installed suites and canonical synthetic-only records. Keep existing PR checkout/jobs/binaries unchanged; commit locally without push.
