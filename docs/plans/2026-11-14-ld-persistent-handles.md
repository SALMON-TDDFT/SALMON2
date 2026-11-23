# LD Persistent Handles Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Retain LD coefficients and per-pole history across calls, preserving public arithmetic, snapshots and atomic failures.

**Architecture:** A private typed native state owns the unchanged shared LD kernel type, coefficients computed once for fixed dt, and rollback buffers. Python keeps only accepted field, totals and clock; pole arrays are exported/imported explicitly for snapshots. step passes one field vector; batch remains a single native call and commits only when all samples are finite. Dynamic IDs, serialized native registry access, explicit close/context manager plus GC cleanup. No fixed machine-wide cap, import-time build, shared-kernel or full-runtime changes.

**Tech Stack:** Fortran ISO_C_BINDING, ctypes/NumPy, existing snapshot contract and audited CMake build.

## Approved design

User approved this design on continuing after the preceding design review. Alternatives: retaining the stateless API leaves allocation/coefficient overhead; Python-only recurrence duplicates the reference kernel and changes the numerical path. Persistent kernel state is selected. Each state has its own dt-dependent coefficients; no model cache or pickling of native resources. LorentzDrude remains data-only and serializable.

## Tasks

1. Add failing contract tests for handles/lifecycle, no pole export on step, old stateless parity (step and batch, multiple pole counts), snapshot rollback and rejected arithmetic. Preserve literal prechange ld.py for benchmark provenance.
2. Add salmonpy/ld_handles.f90: checked allocation of kernel arrays and rollback buffers, scalar finite validation, monotonic int64 IDs, create/advance/export/import/destroy. Keep native/lorentz_drude.f90 byte-identical. Validate every pole and total before commit; batch failure restores starting state.
3. Add salmonpy/_ld_handles.py and adapt ld.py. Finalizer failure cleans up; Python metadata allocated before native mutation; restoration validates against a temporary proxy then imports atomically. Existing stateless symbol retained for comparison. Add CMake packaging and new explicit-library source harness.
4. Audit sources before isolated native build, retain CMake source/object audits, no replacement of existing source binaries. Run focused/full tests, wheel-from-sdist and installed full suite, MPI LD parity.
5. Benchmark independent old/new subprocesses for step/batch and pole counts, setup/time/peak RSS and modeled copies separately, report exact parity and limits. Save records in canonical calc2026, review and commit without push or job changes.

## Completed validation

239 source and installed-wheel tests pass. Source/installed2/3-rank MPI bitwise parity and explicit session cleanup pass; existing4-rank multi-material test passes. All18 old/new benchmark comparisons are bitwise identical. With50000updates and three isolated repetitions, step mode is3.03–3.75× faster and batch is0.999–1.016×. Source-derived pole copies fall168n→48n bytes/call; peak RSS differences are small and do not establish general savings. Shared kernels/runtime patch unchanged, existing binaries preserved, audits retained. Report and artifacts: canonical calc2026/salmonpy/2026-11-14-ld-persistent-handles. Independent review resolved the restore-result allocation gap and found no further blockers.
