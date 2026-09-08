# Wannier90 Fixed-Center Site-Symmetry Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Produce converged symmetry-adapted MLWFs from the global orthonormal LCFO basis without materializing the full affine representation.

**Architecture:** Preserve the 1536-operation affine group as a streamed subspace proof.  Select the maximal full-system fixed-center point subgroup, stream one dense representation matrix at a time into the existing DMN writer, and run Wannier90 with `site_symmetry` enabled.  Reject iteration-limit exhaustion and retain center closure as a publication gate.

**Tech Stack:** Fortran 2008, MPI, EigenExa, spglib, Wannier90 3.1 library mode, SAWF DMN writer, Python source contracts, CMake clean overlays.

---

## Mandatory discipline

For every task, use TDD with a genuine RED, run focused MPI 1/2/4/8 tests,
perform specification and code-quality reviews, resolve every Critical and
Important finding, and build a clean committed-parent overlay.  Do not alter
normal DC LCFO+EigenExa or generalized-eigenvalue Exp-only V3 RT.

### Task SS1: Stream fixed-center representation matrices

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write RED**

Add a small distributed fixed-center group.  Require an API that assembles one
operation into row-owned storage, gathers only that matrix to a designated
writer, matches the existing dense reference, and reports a peak independent
of subgroup order.  Reject invalid operation indices, rank-inconsistent
dimensions, nonfinite input, and allocation overflow.  Forbid a production
`Npoint*Nwann**2` tensor.

**Step 2: Run RED**

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because the one-operation gather API and production stream are
absent.

**Step 3: Implement minimal stream**

Reuse `exchange_dg_point_permuted_orbital_rows` and orbital tiles.  Gather one
row-owned matrix with validated contiguous ownership to rank 0, release it
after the writer append, and retain only measured peak bytes.

**Step 4: Verify and review**

Run construction MPI 1/2/4/8, dense-reference comparisons, route contract,
and `git diff --check`.  Review pullback composition, MPI counts, zero-row
ranks, overflow, and cleanup paths.  Resolve findings and commit SS1 files.

### Task SS2: Publish and validate a fixed-center DMN

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf_dmn.f90` only if a generic streaming defect is found
- Modify: `tests/dg/check_sawf_dmn_format.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write RED**

Require production to call `begin_sawf_dmn`, append every fixed-center
operation exactly once, and call `finish_sawf_dmn` before Wannier90 setup.
For the orthonormal LCFO projection gauge, require `A=I` and
`D_band=D_wann`, with identity, unitarity, covariance, closure, and
process-count-independent fingerprint receipts.  Reject missing inversion and
subgroup order above 48.

**Step 2: Run RED**

```bash
python3 tests/dg/check_sawf_dmn_format.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because the production ground-state route writes no DMN.

**Step 3: Implement minimal DMN publication**

Move fixed-center selection before Wannier90.  Reuse the full atomic catalog
and common inversion center.  Stream each operation matrix to the rank-0
writer as both band and Wannier representations.  Pass identity `A`, current
eigenvalues, and the selected `t_sawf_symop` list.  Atomically finish the DMN;
abort scratch publication on every failure.

**Step 4: Verify and review**

Run DMN format, construction, symmetry projection, and route fixtures on MPI
1/2/4/8.  Review operation ordering, identity-first normalization, affine
translations, off-fragment center, Gamma phases, file cleanup, and workspace
receipts.  Resolve findings and commit SS2 files.

### Task SS3: Enable symmetry-adapted Wannier90 and require convergence

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write RED**

Require `.win` to contain `site_symmetry = .true.` and a strict
`symmetrize_eps`; require a present DMN before library setup.  Add a parser for
the final Wannier90 convergence receipt and reject an exhausted iteration
limit, a still-decreasing final window, missing final state, nonfinite spread,
or absent DMN.  Keep Gamma-real and unitary gates.

**Step 2: Run RED**

```bash
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because setup uses unconstrained mode and validation accepts
the 200-iteration limit.

**Step 3: Implement minimal constrained run**

Write the two symmetry keywords before `wannier_setup`.  Validate the DMN
transaction and parse the rank-0 `.wout` after `wannier_run`; broadcast a
convergence status and iteration count.  Reject rather than falling back to
unconstrained localization.

**Step 4: Verify and review**

Run W90 adapter MPI 1/2/4/8 with real one-band library mode, malformed output
fixtures, route contract, and full-feature incremental build.  Review parser
bounds, stale files, library/file agreement, and error broadcasts.  Resolve
findings and commit SS3 files.

### Task SS4: V3 provenance and genuine Si64 acceptance

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`

**Step 1: Write publication RED**

Require DMN fingerprint/order, site-symmetry flag, convergence iterations,
symmetrization tolerance, fixed-center center/order/inversion, and separate
DMN/W90/full-affine workspace peaks.  Reject missing, zero, nonfinite,
rank-inconsistent, over-limit, unconverged, or non-inversion evidence.

**Step 2: Implement and verify V3**

Extend manifest write/read, expected size, digest, broadcast, replicated
payload validation, rejection code, and production population.  Run checkpoint
MPI 1/2/4/8 and route tests.

**Step 3: Clean overlay and Si64**

Build EigenExa `-j1`, then the full clean overlay with MPI, ScaLAPACK,
EigenExa, spglib, and Wannier90.  Run ideal undisplaced Si64 on 8 MPI ranks and
one OpenMP thread.  Require 384/128 ranks, 1536/32/48 proof orders, converged
symmetry-adapted MLWFs, fixed-center inversion, bounded memory, center closure,
accepted V3, and restart reuse.

**Step 4: Reviews and commit**

Run all affected focused fixtures on MPI 1/2/4/8 and `git diff --check`.
Resolve all Critical/Important findings, perform a clean committed-parent
overlay build, and commit SS4.  Resume AP3 LR/HHG only after this gate passes.
