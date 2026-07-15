# WPW Production API Hardening Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the existing distributed `windowed_kg` operator infrastructure fail closed on lifetime, component-provenance, PP-face, global-candidate, and MPI-failure contract violations.

**Architecture:** Harden the current callback, sparse-builder, and owner-exchange APIs without adding the quadrature assembler or SCF consumer. Preserve the dense solver as a small-system oracle and keep production storage rank-local.

**Tech Stack:** Fortran 2008, MPI, Python 3 contract tests, NumPy fixtures, CMake.

---

### Task 1: Bind callback context safely

**Files:**
- Modify: `src/rt/dg/rt_dg_wpw_matrix_free_adapter.f90`
- Modify: `tests/dg/test_dg_wpw_adapter_mpi.f90`
- Modify: `tests/dg/check_dg_wpw_callback_context.py`

1. Add a failing test requiring a public binder, persistent-target contract,
   bound flag, and rejection of direct unbound callback use.
2. Run the focused test and confirm RED because the binder is absent.
3. Add the minimal binder and make context components private outside the
   module where practical; validate association inputs before setting bound.
4. Bind the MPI fixture through the new API and confirm GREEN.

### Task 2: Separate WW component provenance

**Files:**
- Modify: `src/rt/dg/rt_dg_wpw_matrix_free_adapter.f90`
- Modify: `tests/dg/test_dg_wpw_adapter_mpi.f90`
- Modify: `tests/dg/check_dg_wpw_complete_ww_apply.py`

1. Add a failing test requiring explicit real-local/SIPG and
   complex-nonlocal component tags.
2. Confirm RED because current context accepts untagged arrays.
3. Validate the tags in the binder and H callback; reject ambiguous complete-WW
   plus nonlocal combinations.
4. Confirm the end-to-end H/S dense-equivalence fixture remains GREEN.

### Task 3: Enforce rank-local candidate provenance

**Files:**
- Modify: `src/rt/dg/rt_dg_wpw_sparse_builder.f90`
- Modify: `tests/dg/test_dg_wpw_sparse_builder.f90`
- Modify: `tests/dg/check_dg_wpw_windowed_sparse_builder.py`

1. Add failing fixtures for a PP face candidate and a candidate whose output is
   owned by another rank.
2. Confirm RED because both are currently accepted or filtered silently.
3. Add candidate-origin constants/arrays, require PP volume origin, and reject
   every non-owned input candidate rather than filtering a global list.
4. Confirm the numerical builder fixture and distributed block-action oracle
   are GREEN.

### Task 4: Make MPI collective failure terminal

**Files:**
- Modify: `src/rt/dg/rt_dg_wpw_owner_exchange.f90`
- Modify: `tests/dg/check_dg_wpw_exchange_overflow.py`
- Modify: `tests/dg/test_dg_wpw_owner_exchange_mpi.f90`

1. Add a failing source-contract test requiring one communicator-fatal helper
   after every failed `MPI_Alltoall` and `MPI_Alltoallv`.
2. Confirm RED because current code returns locally.
3. Add a diagnostic helper that calls `MPI_Abort`; retain collective validation
   handshakes only for ordinary rank-local validation failures.
4. Confirm the 3-rank normal and fail-closed ownership fixture is GREEN.

### Task 5: Verify and checkpoint

1. Run all focused Python tests listed in the Task 5 handoff.
2. Compile and run the sparse-builder and weak-form Fortran fixtures.
3. Run the 3-rank owner-exchange and 2-rank adapter fixtures.
4. Run `cmake --build build-mpi-eigenexa -j 1` and `git diff --check`.
5. Request an independent code review and address Critical/Important findings.
6. Commit only this API-hardening checkpoint and confirm `git status --short`
   contains no remaining Task 5 changes.
