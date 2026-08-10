# Symmetry-Adapted LCFO Occupied-Subspace Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the boundary-broken LCFO occupied seed by a memory-bounded, full-system-symmetry-adapted occupied subspace, combine it with the complete buffer `s+p` projector space, and publish only a genuinely closed 384-state Wannier90/V3 route.

**Architecture:** Derive symmetry solely from the full atomic catalog.  Stream fixed-center group actions over distributed fragment-plus-buffer orbital tiles, form a group-averaged occupied projector, select complete invariant clusters totaling rank 128, and take its orthonormal direct sum with the rank-256 complete `s+p` projector space.  Repair the existing affine proof so it measures the actual selected space instead of reporting inferred zeros.

**Tech Stack:** Fortran 2008, MPI, OpenMP, ScaLAPACK/EigenExa, spglib, Wannier90 3.1 library mode, Python source contracts and MPI fixture runners, CMake clean overlays.

---

## Mandatory discipline

For every task, first demonstrate a genuine RED; then run focused MPI 1/2/4/8
verification, specification review, and code-quality review; resolve every
Critical and Important finding; and build a clean-first overlay from the
committed parent plus only that task's prerequisite commits.  Preserve normal
DC LCFO+EigenExa and generalized-eigenvalue Exp-only V3 RT.  Do not loosen
physics tolerances, replace a leaky action by its polar unitary, or fall back
to unconstrained Wannier90.

### Task OA1: Make symmetry closure proof measure the production seed

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write genuine RED**

Add a distributed fixture whose atom/projector permutation is valid but whose
orbital seed has a deliberate component outside its transformed span.  Require
nonzero closure residual, singular values different from one, nonzero measured
workspace, and rejection.  Retain a closed reference that passes identically
on MPI 1/2/4/8.  Add a source contract forbidding proof values assigned from
permutation success or hard-coded zero.

**Step 2: Run RED**

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because the current affine proof reports zero leakage and zero
workspace for the perturbed seed.

**Step 3: Implement streamed residual**

For each operation and orbital tile, apply the real-space pullback to the
actual production seed, assemble `Q^dagger U_g Q`, and measure the orthogonal
residual without storing all transformed orbitals.  Accumulate real allocation
high-water marks.  Use stable reductions and deterministic operation order.

**Step 4: Verify, review, overlay, commit**

Run construction and route fixtures on MPI 1/2/4/8 plus `git diff --check`.
Review adjoint conventions, row ownership, zero-row ranks, allocation overflow,
and cleanup on rejection.  Resolve all Critical/Important findings, perform
the task overlay build, and commit only OA1 files.

### Task OA2: Build the group-averaged occupied projector in bounded memory

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write genuine RED**

For a small fixed-center group and a boundary-perturbed occupied space, require
the streamed averaged projector/eigenproblem to match a dense reference,
reduce closure relative to the input, conserve trace, and keep peak memory
independent of group order.  Test group orders 1 and greater than 1, ranks with
no owned rows, nonfinite data, invalid permutations, overflow, and failure
cleanup.  Forbid production arrays proportional to
`Ngroup*Nocc*Ngrid`.

**Step 2: Run RED**

Run the construction MPI runner and route checker; expect failure because no
occupied-projector averaging API exists.

**Step 3: Implement minimal distributed averaging**

Stream one operation and occupied-orbital tile, exchange required grid rows,
and accumulate the orbit Gram/action matrices.  Solve only the distributed
dense orbit problem.  Return candidate eigenvalues, block metadata, trace,
closure, and measured current/peak bytes without retaining the orbit tensor.

**Step 4: Verify, review, overlay, commit**

Run dense equivalence and adverse fixtures on MPI 1/2/4/8, route checks, and
`git diff --check`.  Review normalization by group order, operation
composition, Gamma reality, MPI reproducibility, and memory accounting.
Resolve all Critical/Important findings, run the clean-first overlay, and
commit OA2.

### Task OA3: Select invariant rank 128 and form the rank-384 seed

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write genuine RED**

Require selection of complete degenerate/invariant blocks totaling the
requested occupied rank and deterministic agreement on MPI 1/2/4/8.  Add RED
cases where rank 128 cuts a cluster, the cluster gap is ambiguous, the
complete `s+p` space loses rank, occupied and projector spaces overlap too
strongly, or the direct sum is not exactly 384.  Require successful cases to
close under every fixed-center operation before DMN writing.

**Step 2: Run RED**

Run construction and route checks; expect failure because production still
uses the raw occupied LCFO states and has no invariant-block gate.

**Step 3: Implement selection and direct sum**

Cluster the averaged-projector spectrum with scale-aware tolerances, verify
each proposed block by group action, and accept exactly 128 states without
splitting a block.  Direct-sum those states with all 256 complete buffer `s+p`
projectors and orthonormalize using the existing generalized algebra.  Measure
rank, conditioning, Gamma-imaginary norm, and fixed-center closure before DMN.

**Step 4: Verify, review, overlay, commit**

Run all success and rejection fixtures on MPI 1/2/4/8, DMN format and route
checks, and `git diff --check`.  Review tolerance scaling, deterministic
tie-breaking, rank logic, preservation of occupied states, and absence of
polar-unitary laundering.  Resolve all Critical/Important findings, run the
clean-first overlay, and commit OA3.

### Task OA4: Publish density, boundary, cluster, and memory evidence

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write genuine RED**

Require V3 fields for original/adapted subspace distance, electron-count
drift, core-interior and buffer-boundary density differences, before/after
closure, selected/rejected eigenvalue edges, cluster gap/block dimensions,
and measured workspace.  Reject missing, zero workspace after nonidentity
work, nonfinite, rank-inconsistent, stale, over-tolerance, or digest-mismatched
receipts.  Test read/write/broadcast/restart on MPI 1/2/4/8.

**Step 2: Implement receipts and gates**

Accumulate density differences on the existing core/buffer masks, reduce them
deterministically, extend the manifest size/digest/broadcast validation, and
connect every rejection to pre-publication cleanup.  Keep interior and
boundary tolerances distinct and named.

**Step 3: Verify, review, overlay, commit**

Run checkpoint, construction, DMN, W90, route, metadata, and Exp-only RT
fixtures on MPI 1/2/4/8 plus `git diff --check`.  Review units, masks,
normalization, restart compatibility, atomic publication, and normal-route
isolation.  Resolve all Critical/Important findings, run the clean-first
full-feature overlay, and commit OA4.

### Task OA5: Genuine Si64 GS, MLWF, and V3 acceptance

**Files:**
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`
- Modify tests or production files only for defects exposed by this acceptance

**Step 1: Establish acceptance RED**

Run the ideal undisplaced Si64 GS with the current committed prerequisite and
capture the existing operation-4 closure rejection as RED.  The run must use
the genuine pseudopotential, 8 MPI ranks, one OpenMP thread, and the clean
full-feature overlay.

**Step 2: Run the adapted production route**

Require normal DC-SCF/LCFO+EigenExa, occupied rank 128, complete `s+p` rank
256, final rank 384, affine orders 1536/32/48, inversion-containing
fixed-center group, strict post-adaptation closure, bounded nonzero workspace,
converged symmetry-adapted Wannier90, center-orbit closure, accepted V3, and
restart reuse.

**Step 3: Review and commit**

Record exact commands, hashes, residuals, density receipts, memory peaks,
Wannier90 iterations/spreads, and restart evidence.  Perform specification and
code-quality reviews, resolve every Critical/Important issue through a fresh
RED task, rerun the clean committed-parent overlay, and commit acceptance
results.

### Task OA6: Resume polarization-primary LR and long-pulse HHG

**Files:**
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/test_si64_harmonic_morphology.py`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_gs.in`
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`

**Step 1: Re-establish RED**

Require the accepted V3 checkpoint and polarization output.  Define the
semilog polarization-spectrum window, harmonic-bin background estimator, and
explicit H2/H4 peak-versus-dip classification before running RT.  Current is
reported only as a secondary consistency observable.

**Step 2: Run LR and HHG**

Use generalized-eigenvalue Exp coefficient RT, the previously justified large
time step, and a sufficiently long pulse to resolve harmonic morphology.  Do
not use displaced atoms or lattice vibration in this acceptance.

**Step 3: Verify, plot, review, commit**

Check inversion selection rules, polarization/current consistency, spectral
resolution, window sensitivity, H2/H4 peak/dip status, even/odd suppression,
and qualitative small-cell limitations.  Produce the requested semilog HHG
figure, run focused MPI 1/2/4/8 fixtures and the final clean-first overlay,
resolve every Critical/Important review finding, update the results document,
and commit.

### Task OA7: Final branch verification and publication

Run every retained DG contract and focused MPI 1/2/4/8 suite, the ideal Si64
GS/restart, LR/HHG analysis, `git diff --check`, and a final clean build from
the committed parent overlay.  Confirm no obsolete DG route or fallback was
reintroduced and normal DC LCFO+EigenExa remains unchanged.  Perform final
specification and code-quality reviews, resolve all Critical/Important
findings, commit any review-only corrections, then push
`codex/wpw-s-orthogonal-complement` to both `origin` and `upstream`.
