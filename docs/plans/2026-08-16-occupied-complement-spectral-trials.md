# Occupied/Complement Spectral Trial Frame Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Preserve the symmetry-adapted occupied LCFO trial block and construct only its complete-s+p-derived orthogonal complement with streamed spectral-basin operators before one full-space Wannier90 call.

**Architecture:** Add a small row-owned trial-frame composition primitive that validates and concatenates preserved and localized complement blocks. Reuse the existing basin preparation, projection, eigensystem, rank-selection, and propagation routines at complement dimension `ntarget-nstate`; change production orchestration rather than duplicating those numerical kernels. Keep one basin operator/eigensystem live at a time and retain the existing single AMN/MMN/DMN and Wannier90 invocation.

**Tech Stack:** Fortran 2008, MPI, EigenExa/LAPACK, Python MPI fixture runners, Wannier90 library interface.

---

### Task 1: Lock the complement-only failure and success contract

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py`

**Step 1: Write the failing test**

Add a `spectral_basin_complement` case with a six-dimensional retained frame,
two preserved occupied directions, and four complement directions. Construct
basin representative candidates that duplicate an occupied direction in the
old full-space scheme. Require the new API to select and propagate exactly four
complement channels and produce a six-dimensional concatenated frame with
Gram-I and occupied/complement cross-Gram defects below tolerance.

Add separate reject cases for a duplicate complement channel and a
rank-disagreeing occupied count.

**Step 2: Run the test to verify it fails**

Run:
`python3 tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py`

Expected: FAIL because the complement composition API does not exist or the
full-space propagation reproduces the duplicate-direction Gram defect.

**Step 3: Commit the RED fixture**

```bash
git add tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py
git commit -m "test(dg): expose occupied spectral channel duplication"
```

### Task 2: Add the row-owned occupied/complement composition primitive

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90`

**Step 1: Implement the minimal API**

Add a public primitive accepting row-owned occupied coefficient rows and
row-owned complement coefficient rows. Collectively agree dimensions and
tolerance, validate exactly-once row ownership, preflight checked allocation
and byte extents, and measure occupied Gram-I, complement Gram-I, and cross-Gram
defects. Return one row-owned concatenated frame and a decomposition-independent
full-projector fingerprint.

Handle `noccupied=0` and `nempty=0` without zero-size reduction/count hazards.

**Step 2: Run the focused test to verify GREEN**

Run:
`python3 tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py`

Expected: PASS for MPI 1/2/4/8, including complement success and reject cases.

**Step 3: Refactor only after GREEN**

Reuse existing checked integer and streamed projector helpers. Do not introduce
a second global retained frame or replicated real-space frame.

**Step 4: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90
git commit -m "feat(dg): compose occupied and complement trial frames"
```

### Task 3: Restrict production basin work to the empty complement

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the failing route assertions**

Require production to derive `nempty=ntarget-nstate`, project prepared basin
operators using only complement state values, select a total rank of `nempty`,
compose occupied and complement rows before AMN/DMN generation, and retain one
`run_dg_w90_gamma_library` call. Reject the obsolete full-`ntarget` basin rank
selection pattern.

**Step 2: Run the route test to verify RED**

Run: `python3 tests/dg/check_dg_overlapping_wannier_route.py`

Expected: FAIL on the current full-space selection and projection calls.

**Step 3: Implement the production orchestration**

- Expose the already orthonormalized complete-s+p-derived complement block from
  the joint retained frame without requesting or materializing DC-LCFO empty
  buffer wavefunctions.
- Prepare/project basin operators at dimension `nempty`.
- Stream one basin orbit at a time and select exactly `nempty` channels.
- Propagate only complement channels.
- Compose the preserved occupied coordinate rows and propagated complement
  rows using the Task 2 primitive.
- Feed the resulting single `spectral_trial_rows` to the existing AMN/DMN and
  one-shot Wannier90 path.
- Remove temporary diagnostics that dump the full selected-rank vector; retain
  concise dimension, Gram defect, workspace, and fingerprint receipts.

**Step 4: Run route checks to verify GREEN**

Run:
`python3 tests/dg/check_dg_overlapping_wannier_route.py`

Run:
`python3 tests/dg/check_obsolete_dg_routes_removed.py`

Expected: PASS.

**Step 5: Commit**

```bash
git add src/gs/main_dft.f90 tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "fix(dg): localize only the empty retained complement"
```

### Task 4: Verify MPI, numerical, memory, and single-Wannier contracts

**Files:**
- Modify only if a genuine regression is found:
  `src/gs/dc/dg_overlapping_wannier_construction.f90`
  `src/gs/main_dft.f90`
  relevant focused tests

**Step 1: Run construction MPI fixtures**

Run:
`python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

Expected: PASS for MPI 1/2/4/8 with invariant fingerprints.

**Step 2: Run EigenExa MPI fixtures**

Run:
`python3 tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py`

Expected: PASS for MPI 1/2/4/8.

**Step 3: Build the overlay target**

Use the existing configured overlay build command recorded by the focused test
runners. Expected: successful compile and link with no new warnings from the
changed routines.

**Step 4: Check source hygiene**

Run: `git diff --check`

Expected: no output.

**Step 5: Commit any test-driven correction separately**

Do not stage user-owned Si64 scripts/input files or the untracked Wannier90
fixture outputs.

### Task 5: Re-run Si64 through the former failure point

**Files:**
- No source changes unless a new RED is first added.

**Step 1: Stop stale calculation processes**

Confirm no previous SALMON, MPI launcher, or Wannier90 process remains before
starting the new binary.

**Step 2: Start the existing Si64 case with unchanged MPI rank count**

Use the corrected executable and the established Si64 input. Capture per-rank
RSS and milestone diagnostics without increasing the MPI rank count or BLAS
thread count.

**Step 3: Verify the former failure boundary**

Require:

- selected complement rank equals `ntarget-nstate`;
- occupied/complement cross-Gram defect is within tolerance;
- complete trial Gram defect is within tolerance;
- only one Wannier90 setup/run appears in the log;
- memory does not contain an all-basin or duplicated spatial trial-frame peak.

**Step 4: Continue through Wannier90 completion**

If it stops, record the exact last milestone, RSS per rank, and diagnostic
receipt. Do not infer memory failure without a measured peak or allocator/OS
evidence.

**Step 5: Commit only source/test changes**

Keep run logs and generated `.wout` files outside the commit.
