# WPW Tail Diagnostic Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the occupied-W outer-shell threshold a non-blocking diagnostic while preserving all structural and physical gates.

**Architecture:** Add a pure tail classifier and retain the existing MPI norm reduction. Separate structural validity from tolerance classification: `total>0` and `0<=outer<=total` (with roundoff allowance) are required collectively, while a valid ratio above tolerance emits a warning and continues to descriptor publication and downstream gates.

**Tech Stack:** Fortran 2008, MPI, Python source-contract tests, existing MPI/EigenExa build.

---

### Task 1: Non-blocking tail classification

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/dg_wpw_wannier_tail_halo.f90`
- Modify: `tests/dg/test_dg_wpw_core_wannier_seed_mpi.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add failing unit tests for `classify_sawf_wannier_buffer_tail`: zero outer,
   below/equal/above tolerance, non-finite values, zero/negative total,
   negative outer, and outer greater than total. Above tolerance must return
   `warning=.true., info=0`; invalid data must return `info/=0`.
2. Run `python3 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`; expect RED
   because the classifier is absent.
3. Implement the pure classifier with a `100*epsilon*max(1,total)` upper-bound
   allowance and no mutation or tolerance change.
4. Add failing source-contract assertions that `buffer_sufficiency` is absent,
   `buffer_tail_validity` synchronizes MPI/classifier failure before descriptor
   construction, and an above-tolerance warning is emitted.
5. Require the existing normalization stages before the warning and the
   metric/Gram, density/electron-count qualification, and fixed-H calls/stages
   after descriptor construction. Verify their ordering relative to the tail
   validity stage.
6. Replace the old tolerance failure with classifier validity plus warning.
   Do not change `dg_wpw_scf_residual_tolerance` or later gates.
7. Run focused tests, full build, and findings-first review.

### Task 2: Fresh continuation evidence

**Files:**
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run a fresh normalized B=6 case past the tail warning.
2. Record the next actual physical gate and its numerical diagnostics.
3. Run final focused tests, full build, diff check, and review.
