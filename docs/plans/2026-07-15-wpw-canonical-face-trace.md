# WPW Canonical-Face Trace Assembly Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Assemble rank-local real-grid WP SIPG face contributions exactly once per canonical face through an explicit trace-provider boundary.

**Architecture:** Add a callback context that supplies both W and P traces at one unwrapped face point, and a scanner that owns only canonical geometry and accumulation. Reuse the existing point evaluator, keep PP face terms structurally absent, and leave remote trace exchange to a later provider/halo checkpoint.

**Tech Stack:** Fortran 2008, existing DG WPW modules, Python source-contract tests, linked Fortran numerical fixtures, CMake, MPI regression fixture.

---

### Task 1: Freeze the trace-provider API

**Files:**
- Create: `src/rt/dg/rt_dg_wpw_face_trace_provider.f90`
- Create: `tests/dg/check_dg_wpw_face_trace_scanner.py`
- Modify: `src/rt/CMakeLists.txt`

**Step 1: Write the failing source-contract test**

Require a provider context with an associated callback, explicit bind/unbind
operations, face identity and unwrapped grid arguments, both-side W/P values
and gradients, and a status result. Reject `mpi_`, dense H/S names, and PP face
outputs.

**Step 2: Run the test to verify RED**

Run: `python3 tests/dg/check_dg_wpw_face_trace_scanner.py`

Expected: FAIL because `rt_dg_wpw_face_trace_provider.f90` is absent.

**Step 3: Implement the minimal provider context**

Use an abstract interface with a caller-owned `class(*)` context pointer and
procedure pointer. Validate association in a single `evaluate_wpw_face_traces`
entry point. Binding must not copy or retain a procedure-local target.

**Step 4: Run the test to verify GREEN**

Run: `python3 tests/dg/check_dg_wpw_face_trace_scanner.py`

Expected: PASS provider API contract.

**Step 5: Commit**

```bash
git add src/rt/CMakeLists.txt src/rt/dg/rt_dg_wpw_face_trace_provider.f90 tests/dg/check_dg_wpw_face_trace_scanner.py
git commit -m "feat: add WPW face trace provider API"
```

### Task 2: Add canonical-face geometry scanning

**Files:**
- Create: `src/rt/dg/rt_dg_wpw_face_trace_scanner.f90`
- Modify: `tests/dg/check_dg_wpw_face_trace_scanner.py`
- Modify: `src/rt/CMakeLists.txt`

**Step 1: Extend the failing contract test**

Require `assemble_wpw_canonical_face_grid`, canonical `(K-,K+,axis,side)`
validation, `normal(axis)=side`, `h_normal=hgs(axis)`, tangential area weight,
and a call to `assemble_wpw_canonical_face_point`. Reject loops over
`1:n_frag`, hidden MPI, PP face storage, and dense matrices.

**Step 2: Run the test to verify RED**

Run: `python3 tests/dg/check_dg_wpw_face_trace_scanner.py`

Expected: FAIL because the scanner entry point is absent.

**Step 3: Implement minimal geometry and accumulation**

Validate `K-<K+`, axis 1..3, side ±1, positive spacing, sorted ids, output
shape, and intersecting face extents. Loop over only the two tangential grid
indices. At each point call the bound provider and then
`assemble_wpw_canonical_face_point`; accumulate into a temporary block and
publish it only after the entire face succeeds.

**Step 4: Run the test to verify GREEN**

Run: `python3 tests/dg/check_dg_wpw_face_trace_scanner.py`

Expected: PASS scanner source contract.

**Step 5: Commit**

```bash
git add src/rt/CMakeLists.txt src/rt/dg/rt_dg_wpw_face_trace_scanner.f90 tests/dg/check_dg_wpw_face_trace_scanner.py
git commit -m "feat: scan canonical WPW face quadrature"
```

### Task 3: Prove numerical face accumulation

**Files:**
- Create: `tests/dg/test_dg_wpw_face_trace_scanner.f90`
- Create: `tests/dg/run_dg_wpw_face_trace_fixture.py`
- Modify: `tests/dg/check_dg_wpw_face_trace_scanner.py`

**Step 1: Write the failing linked fixture**

Construct a two-fragment periodic mock with one face normal to x and a 2x2
tangential face. Bind a deterministic provider whose traces depend on the
unwrapped grid point. Compute the expected integral by four direct calls to
the existing point evaluator and assert the scanner result, point count,
normal, face identity, and grid coordinates.

**Step 2: Run the fixture to verify RED**

Run: `python3 tests/dg/run_dg_wpw_face_trace_fixture.py --build-dir build-mpi-eigenexa`

Expected: FAIL to compile or link because scanner behavior is incomplete.

**Step 3: Complete only the behavior required by the fixture**

Fix geometry indexing and provider output validation without adding remote
communication or PP face paths.

**Step 4: Run the fixture to verify GREEN**

Run: `cmake --build build-mpi-eigenexa -j2`

Run: `python3 tests/dg/run_dg_wpw_face_trace_fixture.py --build-dir build-mpi-eigenexa`

Expected: `[100%] Built target salmon` and `PASS canonical WPW face trace numerical fixture`.

**Step 5: Commit**

```bash
git add tests/dg/test_dg_wpw_face_trace_scanner.f90 tests/dg/run_dg_wpw_face_trace_fixture.py tests/dg/check_dg_wpw_face_trace_scanner.py src/rt/dg/rt_dg_wpw_face_trace_scanner.f90
git commit -m "test: verify canonical WPW face trace quadrature"
```

### Task 4: Add fail-closed and periodic uniqueness cases

**Files:**
- Modify: `tests/dg/test_dg_wpw_face_trace_scanner.f90`
- Modify: `src/rt/dg/rt_dg_wpw_face_trace_scanner.f90`

**Step 1: Add failing numerical cases one at a time**

Cover provider failure after a partial face, NaN trace, P-trace mismatch,
invalid axis/side/face identity, empty tangential intersection, and both
distinct wrapped faces in the two-fragment periodic case. Assert a zero output
for every failure.

**Step 2: Run each new case to verify RED**

Run after each addition:
`python3 tests/dg/run_dg_wpw_face_trace_fixture.py --build-dir build-mpi-eigenexa`

Expected: FAIL at the newly added assertion.

**Step 3: Implement the minimal validation**

Accumulate into a temporary matrix, check finite callback outputs before the
point evaluator, and assign the caller output only after a successful complete
scan.

**Step 4: Run the fixture to verify GREEN**

Expected: all numerical and failure cases PASS.

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_wpw_face_trace_scanner.f90 tests/dg/test_dg_wpw_face_trace_scanner.f90
git commit -m "fix: harden WPW face trace assembly"
```

### Task 5: Focused regression and review checkpoint

**Files:**
- Modify only files required by review findings.

**Step 1: Run focused source contracts**

Run: `for f in tests/dg/check_dg_wpw_*.py; do python3 "$f" || exit 1; done`

Expected: every WPW contract prints PASS.

**Step 2: Run build and numerical fixtures**

Run: `cmake --build build-mpi-eigenexa -j2`

Run: `python3 tests/dg/run_dg_wpw_grid_point_fixture.py --build-dir build-mpi-eigenexa`

Run: `python3 tests/dg/run_dg_wpw_face_trace_fixture.py --build-dir build-mpi-eigenexa`

Expected: all PASS.

**Step 3: Run the distributed regression fixture**

Run: `python3 tests/dg/run_dg_wpw_adapter_mpi_fixture.py --build-dir build-mpi-eigenexa`

Expected: MPI equivalence PASS with a near-roundoff error.

**Step 4: Request findings-first code review**

Review canonical uniqueness, unwrapped periodic coordinates, callback lifetime,
all-or-zero failure behavior, absence of PP face terms, and absence of global
scans/storage. Apply each valid finding with a new RED test first.

**Step 5: Verify and commit review fixes**

Run: `git diff --check`

Run all commands from Steps 1-3 again, then commit only scoped changes.
