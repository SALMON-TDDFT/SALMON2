# Si8 MPI8 Memory Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add and run a reproducible eight-rank Si8 overlapping-Wannier ground-state memory diagnostic without asynchronous process suspension.

**Architecture:** Materialize a tracked Si8 conventional-cell fixture from existing silicon data, then run the committed clean-overlay binary through a monitored MPI8 wrapper.  The wrapper uses unbuffered output and periodic read-only process/memory sampling, with one-shot whole-job termination at safety thresholds.

**Tech Stack:** SALMON Fortran 2008, Open MPI, Wannier90, EigenExa, Python 3, macOS `ps`, `vm_stat`, and `sysctl`.

---

### Task 1: Add the Si8 fixture contract

**Files:**
- Create: `tests/dg/data/si8_overlapping_wannier/inputfile.in`
- Create: `tests/dg/data/si8_overlapping_wannier/atom.dat`
- Create: `tests/dg/check_si8_overlapping_wannier_fixture.py`

**Step 1: Write the failing fixture checker**

Require actual silicon (`izatom=14`), eight atoms, 32 electrons, a `16^3`
grid, `2 x 2 x 2` fragments/rank-grid decomposition, overlapping-Wannier mode,
and no unresolved placeholders.

**Step 2: Run the checker and verify RED**

Run: `python3 tests/dg/check_si8_overlapping_wannier_fixture.py`

Expected: FAIL because the fixture files do not exist.

**Step 3: Materialize the minimal fixture**

Use the conventional-cell coordinates from
`samples/benchmark/bulk_Si/input_sc_Si.dat` and the reviewed overlapping-Wannier
controls from the larger fixture.  Reference `Si_rps.dat` in the run directory.

**Step 4: Run the checker and verify GREEN**

Run: `python3 tests/dg/check_si8_overlapping_wannier_fixture.py`

Expected: PASS.

### Task 2: Add a nonintrusive MPI memory monitor

**Files:**
- Create: `tests/dg/run_si8_overlapping_wannier_memory.py`
- Test: `tests/dg/test_si8_memory_monitor.py`

**Step 1: Write RED tests**

Test parsing of `vm_stat` and `ps`, phase extraction from a growing log,
per-rank and available-memory thresholds, and whole-process-group termination.
Assert that the implementation contains no `SIGSTOP` or `SIGCONT`.

**Step 2: Run tests and verify RED**

Run: `python3 tests/dg/test_si8_memory_monitor.py`

Expected: FAIL because the runner is absent.

**Step 3: Implement the runner**

Create a fresh persistent run directory, copy the fixture and Si
pseudopotential, record SHA-256 provenance, launch exactly eight MPI ranks with
OMP/BLAS threads set to one and unbuffered Fortran output, sample every 30
seconds, and terminate the MPI process group once on a configured safety gate.

**Step 4: Run tests and verify GREEN**

Run:

```bash
python3 tests/dg/test_si8_memory_monitor.py
python3 tests/dg/check_si8_overlapping_wannier_fixture.py
git diff --check
```

Expected: PASS.

### Task 3: Run the committed clean binary on MPI8

**Files:**
- Use: `/Users/otobetoshihito/SALMON-dev/verification/20260813-generator-proof-build-boz/salmon`
- Create: `/Users/otobetoshihito/SALMON-dev/verification/20260813-si8-memory-mpi8/`

**Step 1: Verify the binary and configuration**

Record executable hash and confirm MPI, EigenExa, spglib, and Wannier90 are
enabled in the clean build cache.

**Step 2: Run the monitored MPI8 calculation**

Run the new runner with an 8 GiB available-memory floor and a conservative
per-rank RSS ceiling.  Do not run any other production calculation concurrently.

**Step 3: Evaluate evidence**

Require the phase log to reach Wannier90 and the post-gauge proof or identify
the exact preceding marker.  Plot or summarize total/per-rank RSS versus time
and distinguish a bounded peak from monotonic growth.

**Step 4: Decide the Si64/C64 next step**

Only if Si8 completes with bounded memory should the larger case be redesigned
or resumed.  If Si8 grows or stops, diagnose that exact phase before any rerun.
