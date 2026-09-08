# Wannier90 Convergence Tolerance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the generated Wannier90 convergence tolerance numerically attainable for the Si64 spread scale without weakening SALMON's explicit-convergence validation.

**Architecture:** Keep Wannier90 as the sole convergence authority. Change only the generated `conv_tol` from `1.d-12` to `1.d-10`, preserve the five-step window and 200-iteration limit, and lock the generated input contract with a focused route test.

**Tech Stack:** Fortran 2008, Wannier90 library interface, MPI fixture runner, Python route checks.

---

### Task 1: Add the generated-input regression

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the failing test**

Require the production source to contain `conv_tol = 1.d-10` and reject the obsolete `conv_tol = 1.d-12` literal in the Gamma-library setup path.

**Step 2: Run test to verify it fails**

Run: `python3 tests/dg/check_dg_overlapping_wannier_route.py`

Expected: FAIL because the generated input still contains `conv_tol = 1.d-12`.

**Step 3: Commit the RED**

```bash
git add tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "test(dg): require attainable Wannier90 tolerance"
```

### Task 2: Change the generated tolerance

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Test: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Implement the minimal change**

Change the generated line to:

```fortran
write(unit,'(a)')'conv_tol = 1.d-10'
```

Keep `conv_window = 5`, `num_iter = 200`, and `validate_dg_w90_convergence_log` unchanged.

**Step 2: Run the focused test**

Run: `python3 tests/dg/check_dg_overlapping_wannier_route.py`

Expected: PASS.

**Step 3: Run focused MPI verification**

Run the existing Wannier90 MPI fixture for 1, 2, 4, and 8 ranks from a temporary directory, then run the obsolete-route check and `git diff --check`.

Expected: all PASS.

**Step 4: Build production**

Run: `cmake --build /tmp/salmon-wpw-full --target salmon -j1`

Expected: build completes successfully.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_w90.f90
git commit -m "fix(dg): use attainable Wannier90 convergence tolerance"
```

### Task 3: Re-run Si64 through Wannier90

**Files:**
- Inspect: `/Users/otobetoshihito/SALMON-dev/verification/20260816-si64-direct-retained-frame-rerun-mpi8`
- Create: a new timestamped verification directory outside the worktree

**Step 1: Verify no prior SALMON/MPI calculation remains active**

Use a read-only process listing and do not terminate unrelated processes.

**Step 2: Launch the same Si64 MPI-8 input with the rebuilt executable**

Keep MPI rank count and thread settings identical to the prior run.

**Step 3: Monitor bounded checkpoints**

Confirm DC-SCF, DC-LCFO, direct retained frame, translation sectors, Wannier90 convergence, and the first post-Wannier stage. Record RSS per rank and convergence iteration.

**Step 4: Report evidence**

Report whether the explicit Wannier90 convergence sentence appears before iteration 200 and whether SALMON passes the former validation failure. Do not claim full success unless the run and focused tests provide fresh evidence.

