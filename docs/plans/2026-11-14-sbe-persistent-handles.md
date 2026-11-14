# Persistent SBE handles implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove repeated Python/Fortran full-model and density transfers from SBE stepping while preserving the existing SBE API, numerical update order, independent states and snapshots.

**Architecture:** Add a typed internal handle module in salmonpy/sbe_handles.f90 using unchanged bundled sbe_types/velocity kernels. One immutable native model copy per SBEInitialState is cached lazily; states hold independent solver/trial density buffers and accepted A/J. A step copies rho only into the native trial buffer, runs unchanged kernels, validates every scalar density component/current, then swaps buffers atomically. Python sends A3/receives J3 and exports rho only on density/snapshot access. Explicit close and weakref finalizers release state handles; model helper outlives all states. All native registry calls are serialized. Existing stateless ABI remains for a literal legacy baseline. No full runtime patch/shared-kernel edits.

**Tech Stack:** Fortran ISO_C_BINDING, ctypes/NumPy, existing snapshot validation, CMake kernel build with retained source/object audits.

## Decisions and limits

The user continued after public ABI1 and accepted independent native handles. This advances persistent SBE implementation (list19), not native TDDFT reentrancy. A shared native model avoids copying p/rvnl/eigen/weights and rebuilding delta/occup every call. Two rho buffers preserve failure atomicity; one native rho shadow copy per step remains. Resident memory may rise for many idle states, so measure rather than claim unconditional reduction. No arbitrary machine-wide memory/state limit; sizes checked for integer overflow and allocations return errors. Scalar IEEE inquiries only. Outer-k kernel OpenMP loops/arithmetic remain byte-identical.

## Task1: Red contract tests

Create test_sbe_handles.py covering shared native model, vector-only steps, lazy immutable density export, independent/closed state lifecycle, weakref cleanup, exact legacy-current/density parity across4orders×2vnl settings, snapshot/restore and rejected numerical steps. Preserve the prechange sbe.py for benchmark provenance. Do not overwrite the previously validated source native binary.

## Task2: Native/Python implementation

Add handle module with dynamic linked registries, positive monotonic IDs, model/state ownership, status returns, scalar finite loops, checked dimensions/allocations and transactional rho swap. Cache model helper externally in a weak-key map; retain data-only initial-model serialization and frozen-array reconstruction. Add Python adapter/finalizers; restore validates token into a temporary proxy, imports native rho atomically, then updates Python counters/A/J. Retain exact public charge/matter signs.

## Task3: Audited isolated build and correctness

Build/install new kernel library in a new /tmp prefix with existing CMake prohibition targets. Source tests explicitly inject this verified new library via a validation harness instead of replacing source/_native. Run full tests/golden fixtures and literal stateless baseline from the same newly audited library. Existing unchanged executables are not rebuilt/replaced. No GS collection.

## Task4: Performance/distribution

Independent subprocess benchmarks compare setup/step time, modeled explicit copy bytes and peak RSS for several band/k/state counts. Numerical differences must be reported separately from timings. Update guides and package source/module. Build sdist→wheel and run installed suite outside checkout plus MPI point-distribution smoke. Review, canonical calc2026 records, commit without push or production-job changes.
