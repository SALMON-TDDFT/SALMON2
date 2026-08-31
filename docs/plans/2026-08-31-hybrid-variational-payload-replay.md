# Hybrid Variational Payload Replay Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Capture the eight-rank Hybrid fixed-payload inputs once and replay their extent validation without repeating Si64 DC and Wannier construction.

**Architecture:** Add a versioned, rank-sharded diagnostic bundle beside the existing variational-payload validation. Capture is opt-in through `SALMON_DG_VARIATIONAL_PAYLOAD_CAPTURE`; a standalone MPI test executable reads the shards, reports exact observed extents, and invokes the unchanged production freeze routine.

**Tech Stack:** Fortran 2008 stream I/O, MPI, Python test runner, CMake SALMON build.

---

### Task 1: Specify and implement payload bundle round trip

**Files:**
- Modify: `src/gs/dc/dg_hybrid_variational_payload.f90`
- Create: `tests/dg/test_dg_hybrid_variational_payload_replay_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_variational_payload_replay_mpi.py`

**Step 1: Write the failing MPI round-trip test**

Construct the existing four-row distributed fixture.  Call new public
`write_dg_hybrid_variational_payload_bundle` and
`read_dg_hybrid_variational_payload_bundle` entry points.  Require exact
preservation of global count, row IDs, all four matrices, and all three
fingerprints on 1, 2, and 4 ranks.  Require an existing target prefix to be
rejected rather than overwritten.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/run_dg_hybrid_variational_payload_replay_mpi.py
```

Expected: compilation fails because the writer and reader do not exist.

**Step 3: Implement the minimal bundle format**

In `dg_hybrid_variational_payload` add a format version constant and the two
public routines.  Each shard is an unformatted stream containing:

```text
version, rank_count, rank, global_basis_count, local_row_count,
basis_fingerprint, metric_fingerprint, interface_fingerprint,
row_ids, metric_rows, kinetic_rows, nonlocal_rows, interface_rows
```

Write `<prefix>.rankNNNNNN.tmp`, close it, and publish it with the same
same-directory `rename` pattern used by
`dg_overlapping_wannier_checkpoint.f90`.  After a barrier, rank zero publishes
`<prefix>.manifest`.  Use `status='new'` and reject existing final paths.
The reader validates version, rank count, rank identity, nonnegative allocation
sizes, complete reads, and collective metadata agreement before allocating the
arrays.

**Step 4: Run GREEN**

Run the replay MPI runner.  Expected: PASS on 1, 2, and 4 ranks.

**Step 5: Commit**

Stage only these bundle and test files and commit:

```text
test(dg): add variational payload replay bundle
```

### Task 2: Capture the failing Si64 boundary

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`

**Step 1: Write RED route assertions**

Require `main_dft` to read `SALMON_DG_VARIATIONAL_PAYLOAD_CAPTURE` and call the
bundle writer immediately before `freeze_dg_hybrid_variational_payload` only
when the value is nonempty.  Require capture failure to stop before payload
freeze so a partial diagnostic cannot be mistaken for evidence.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/check_dg_hybrid_continuation_route.py
```

Expected: FAIL because no capture hook exists.

**Step 3: Implement the opt-in hook**

Import the writer, resolve the environment variable once in the continuation
branch, and pass the exact arguments already passed to freeze.  Emit on rank
zero:

```text
[HYBRID-VARIATIONAL-PAYLOAD-CAPTURE] prefix=<path>
```

Do nothing when the environment variable is unset.

**Step 4: Run GREEN and build**

Run:

```text
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/run_dg_hybrid_variational_payload_replay_mpi.py
python3 tests/dg/run_dg_hybrid_variational_payload_mpi.py
cmake --build build-hybrid-commit -j2
```

Expected: PASS and successful SALMON build.

**Step 5: Commit**

Stage only capture-hook hunks and commit:

```text
test(dg): capture Hybrid variational payload
```

### Task 3: Capture once, replay, and identify the extent defect

**Files:**
- Preserve generated evidence under `verification-si64-dg-continuation-20260830/` without staging it.
- Modify production code only after replay proves the defect.

**Step 1: Run one fresh capture calculation**

Use a new output directory and no timeout:

```text
SALMON_DG_VARIATIONAL_PAYLOAD_CAPTURE=<absolute-prefix> \
OMP_NUM_THREADS=1 python3 tests/dg/run_dg_hybrid_si64_continuation_rt.py \
  build-hybrid-commit/salmon \
  verification-si64-dg-continuation-20260830/si64-run-payload-capture-20260831 \
  --ranks 8
```

Preserve the GS log and all eight shards whether the run passes or fails.

**Step 2: Replay without DC**

Run the replay executable on eight ranks against the captured prefix.  Print
for every rank:

```text
rank, global_basis_count, local_row_count, row_id_min, row_id_max,
metric_shape, kinetic_shape, nonlocal_shape, interface_shape
```

Expected: reproduce `invalid variational fixed payload extent` and identify the
specific failed predicate.

**Step 3: Add the minimal failing regression**

Reduce the captured condition to the smallest MPI fixture in the appropriate
focused test.  Run it RED before changing production code.

**Step 4: Fix only the proven defect**

Make the smallest source correction, rerun the replay and affected focused
tests, and commit the fix separately.

**Step 5: Resume the parent handoff plan**

Rerun all protected verification and the final fresh eight-rank Si64 GS-to-RT
calculation required by
`docs/plans/2026-08-30-hybrid-dc-symmetry-handoff.md`.  Do not treat replay
success as final Si64 acceptance.
