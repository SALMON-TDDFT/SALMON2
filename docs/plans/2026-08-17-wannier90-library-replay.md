# Wannier90 Library Replay Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Export the assembled Gamma-point Wannier90 library inputs as a standard replay bundle and run that bundle with standalone `wannier90.x` without repeating SALMON preprocessing.

**Architecture:** Add a small MPI-safe export helper around the existing standard `.eig/.amn/.mmn` writer and call it immediately before the library solve when an environment variable requests export. Add a Python replay driver that works in an isolated temporary directory and preserves the exported bundle. Keep the normal production path byte-for-byte inactive when export is disabled.

**Tech Stack:** Fortran 2008, MPI, Wannier90 standard text formats, Python 3 `unittest`/`subprocess`, existing SALMON focused MPI test harnesses.

---

### Task 1: Lock the standard Gamma payload ordering with a failing test

**Files:**
- Modify: `tests/dg/test_sawf_local_seed_writer.py`

**Step 1: Write the failing test**

Extend the generated Fortran fixture with a two-band, two-projection,
two-neighbor payload whose complex values are unique by index.  Parse the
resulting `.amn` and `.mmn` records in Python and assert the exact Wannier90
ordering, including reciprocal vectors from `neighbor_gvec`.

**Step 2: Run the test to verify current compatibility**

Run: `python3 tests/dg/test_sawf_local_seed_writer.py`

Expected: PASS if the existing writer is directly reusable.  If it fails, the
failure must identify an ordering mismatch before production integration.

**Step 3: Commit**

```bash
git add tests/dg/test_sawf_local_seed_writer.py
git commit -m "test: lock Wannier90 replay matrix ordering"
```

### Task 2: Add an MPI-safe replay export helper

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_w90_mpi.py`

**Step 1: Write failing MPI tests**

Add fixture modes covering:

- disabled export leaves no files;
- enabled export writes `.eig/.amn/.mmn` with exact payload;
- `.win/.dmn` are copied into the bundle;
- missing directory and nonfinite data reject collectively;
- rank-disagreeing directory/enable state rejects without hanging;
- bundle fingerprints match on MPI 1/2/4/8.

Run: `python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py`

Expected: FAIL because the export helper does not exist.

**Step 2: Implement the helper**

Add a public routine with explicit communicator, directory, seed, source
`.win/.dmn`, root-owned eigen/AMN/MMN arrays, and neighbor vectors.  It must:

1. agree enabled state and directory bytes before branching;
2. validate checked dimensions and finite payload collectively;
3. call `write_sawf_local_eig_amn_mmn` only on rank zero;
4. copy `.win/.dmn` using checked Fortran stream I/O;
5. verify all five outputs are non-empty;
6. broadcast root status/message and clean partial outputs on failure.

No dense payload copy or all-rank gather is allowed.

**Step 3: Run the focused MPI test**

Run: `python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py`

Expected: PASS on MPI 1/2/4/8.

**Step 4: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_w90.f90 \
  tests/dg/test_dg_overlapping_wannier_w90_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_w90_mpi.py
git commit -m "feat: export Wannier90 library replay bundle"
```

### Task 3: Integrate opt-in export before the production library call

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Add failing route assertions**

Assert that production reads `SALMON_DG_W90_REPLAY_DIRECTORY`, calls the export
helper after `assemble_dg_w90_gamma_matrices` and before
`run_dg_w90_gamma_library`, passes `w90_nncell`, and leaves the library call
unconditional.

Run: `python3 tests/dg/check_dg_overlapping_wannier_route.py`

Expected: FAIL because the hook is absent.

**Step 2: Implement production integration**

Read the environment variable into a fixed-length buffer with checked status,
collectively agree its exact trimmed value, and invoke the helper only when it
is non-empty.  Use the current seed and setup-produced `.win/.dmn` paths.  On
failure, use the existing collective fatal diagnostic path.  Do not retain any
new arrays after export.

**Step 3: Verify route and build**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
cmake --build /tmp/salmon-wpw-full --target salmon -j2
```

Expected: route PASS and SALMON links successfully.

**Step 4: Commit**

```bash
git add src/gs/main_dft.f90 tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "feat: expose opt-in Wannier90 replay export"
```

### Task 4: Add the standalone replay driver

**Files:**
- Create: `tests/dg/replay_dg_wannier90_bundle.py`
- Create: `tests/dg/test_replay_dg_wannier90_bundle.py`

**Step 1: Write failing Python tests**

Test required-file validation, immutable source handling, `.win`
`symmetrize_eps` replacement/insertion, executable failure propagation, and
successful capture of `.wout` using a fake executable.

Run: `python3 tests/dg/test_replay_dg_wannier90_bundle.py`

Expected: FAIL because the driver is absent.

**Step 2: Implement the minimal driver**

Use `argparse`, `tempfile`, `shutil`, and `subprocess.run` without shell
evaluation.  Accept bundle directory, seed, executable, output directory, and
optional `--symmetrize-eps`.  Copy required files into a temporary directory,
edit only the copy, run `wannier90.x <seed>`, then copy `.wout` and a small JSON
receipt containing command, return code, elapsed time, and input hashes to the
requested output directory.

**Step 3: Run tests**

Run: `python3 tests/dg/test_replay_dg_wannier90_bundle.py`

Expected: PASS.

**Step 4: Commit**

```bash
git add tests/dg/replay_dg_wannier90_bundle.py \
  tests/dg/test_replay_dg_wannier90_bundle.py
git commit -m "test: add standalone Wannier90 replay driver"
```

### Task 5: End-to-end smoke and handoff

**Files:**
- Modify: `docs/plans/2026-08-17-wannier90-library-replay-design.md` only if observed behavior requires clarification

**Step 1: Run focused verification**

```bash
python3 tests/dg/test_sawf_local_seed_writer.py
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/test_replay_dg_wannier90_bundle.py
git diff --check
```

Expected: all PASS and no whitespace errors.

**Step 2: Export one real Si64 bundle**

Create an empty replay directory, run Si64 with MPI 8, `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, and `SALMON_DG_W90_REPLAY_DIRECTORY` set.  Stop after
the bundle is confirmed complete if the in-process solve remains long.

Expected: all five replay inputs exist and memory does not increase by a dense
matrix copy.

**Step 3: Replay standalone**

Run the replay driver with the bundled patched `wannier90.x` and retain its
`.wout`.  Compare site-symmetry residual behavior without rerunning SALMON.

**Step 4: Final review and commit any documentation correction**

```bash
git add docs/plans/2026-08-17-wannier90-library-replay-design.md
git commit -m "docs: record Wannier90 replay verification"
```
