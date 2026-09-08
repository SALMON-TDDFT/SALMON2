# DG Continuation Coding-Rules Compliance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the Task 1--5 DG continuation code comply with the SALMON Fortran coding rules without changing numerical behavior.

**Architecture:** Add one focused static checker, observe it fail on the current code, then apply mechanical naming, import, and line-length changes only.  Preserve direct MPI where no equivalent SALMON wrapper satisfies the datatype and recoverable-error requirements.

**Tech Stack:** Fortran 2008, MPI, Python 3, CMake

---

### Task 1: Add the coding-rules regression

**Files:**
- Create: `tests/dg/check_dg_continuation_coding_rules.py`

1. Check the Task 1--5 Fortran files for lines longer than 132 columns.
2. Reject unrestricted `use mpi` imports.
3. Reject new communicator/rank declarations named `comm` or `rank`.
4. Run the checker and confirm that it fails on the current implementation.

### Task 2: Apply the minimum compliance refactor

**Files:**
- Modify: `src/common/dg_hybrid_continuation_residuals.f90`
- Modify: `src/common/dg_hybrid_sparse_operators.f90`
- Modify: `src/gs/dc/dg_hybrid_continuation_controller.f90`
- Modify: `src/gs/dc/dg_hybrid_continuation_scf.f90`
- Modify: `src/gs/dc/dg_hybrid_continuation_state.f90`
- Modify: `src/gs/dc/dg_hybrid_sipg_operator.f90`
- Modify the five corresponding Fortran MPI fixtures under `tests/dg/`.

1. Rename only parallel identifiers; do not change procedure ordering or data flow.
2. Replace each unrestricted MPI import with `use mpi, only: ...`.
3. Wrap the single overlength fixture line.
4. Run the static checker and confirm PASS.

### Task 3: Verify and commit

1. Run the Task 1--5 Python MPI runners.
2. Compile the changed production modules in a no-MPI configuration.
3. Stage only the two plan files, checker, and Task 1--5 compliance files.
4. Run `git diff --cached --check` and inspect `git diff --cached`.
5. Commit as `style(dg): follow SALMON coding rules`.
