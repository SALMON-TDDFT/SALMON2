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

Select a deterministic generating set from the complete affine product table.
For each generator and orbital tile, apply the real-space pullback to the
actual production seed, assemble `Q^dagger U_g Q`, and measure the orthogonal
residual without storing all transformed orbitals.  Assemble the EigenExa
metric in row tiles using its padded local dimensions, never scalar
collectives.  Accumulate real allocation high-water marks.  Use stable
reductions and deterministic operation order.  Retain all-operation atomic
and cocycle provenance because generator invariance proves closure of the
generated full affine group.

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

### Task OA5a: Reproduce and lock the core-first affine failure

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `src/gs/main_dft.f90` only for bounded diagnostic publication

**Step 1: Write the failing regression**

Construct two translated fragments with overlapping buffers.  Put one smooth
orbital and one complete atomic projector across their shared face.  Assert
that the existing core-first restriction gives a nonzero translation residual
while partition composition gives the exact reference function.  Require
separate occupied, projector, and direct-sum residual labels so the genuine
Si64 failure cannot again be misattributed to only the occupied block.

**Step 2: Verify genuine RED**

Run `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py` on MPI
1/2/4/8 and `python3 tests/dg/check_dg_overlapping_wannier_route.py`.  Expect
the new composition assertion to fail because no buffer-first composer exists.
Retain the fresh ideal-Si64 evidence values `4.93771`, `8.10362`, and
`10.0621` in the design/results record, not as fixture golden tolerances.

**Step 3: Review, overlay, and commit**

Review that the RED fails for fragment-face truncation rather than MPI order,
floating-point noise, or a relaxed threshold.  Resolve every
Critical/Important test finding, run `git diff --check`, perform a clean-first
full-feature build from committed parent `5c1d0f5d`, and commit the regression.

### Task OA5b: Stream buffer-composed physical-grid orbital tiles

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Extend RED contracts**

Require a composer whose inputs are fragment-plus-buffer physical IDs,
partition weights, and an orbital tile.  It must route equal physical IDs to
a deterministic owner, sum `partition_weight*orbital_value`, reject
missing/excess partition coverage, duplicate output ownership, nonfinite
values, rank-inconsistent tile shapes, and integer/count overflow.  Require a
measured nonzero workspace receipt and rank-independent output fingerprint.

**Step 2: Implement the minimal streamed composer**

Add `compose_dg_buffered_orbital_tile_to_physical_grid`.  Use bounded orbital
tiles and MPI count/displacement helpers; do not broadcast a complete orbital,
replicate the global grid, or retain all fragment buffers.  Return sorted
owned physical IDs and values so the existing point-action exchange can act
without fragment semantics.  Keep normal DC LCFO+EigenExa untouched.

**Step 3: Focused verification**

Run the construction fixture on MPI 1/2/4/8, including reversed rank order,
uneven physical-grid ownership, a translation crossing the fragment face,
and adverse partition/shape/nonfinite cases.  Run the route and obsolete-route
contracts and `git diff --check`.

**Step 4: Specification review, quality review, overlay, commit**

Review physical normalization (`sum_f w_f psi_f` for composition, while
`sqrt(w_f)` is only the equivalent replicated-slab inner-product form),
deterministic ownership, collective error agreement, 64-bit counts, memory scaling, and
cleanup on every failure.  Resolve every Critical/Important finding with a
new RED.  Build the committed-parent prerequisite with SPGLIB, EigenExa,
ScaLAPACK, MPI, and Wannier90 using `cmake --build <overlay> --clean-first -j1`,
rerun MPI 1/2/4/8 focused verification, and commit OA5b.

### Task OA5c: Move symmetry adaptation and affine proof onto composed seeds

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write production-order RED**

Require LCFO and complete `s+p` values to remain on fragment-plus-buffer slabs
until composed.  Assert that composition precedes fixed-center averaging,
rank-384 orthonormalization, full-affine generator residuals, DMN creation, and
Wannier90.  Forbid `accumulate_*_to_core` or `core_mask` as the support of any
pre-Wannier symmetry proof.

**Step 2: Integrate buffer-first data flow**

Compose occupied and projector orbital tiles to distributed physical-grid
owners, then run the existing fixed-center group average and complete-cluster
selection on that representation.  Form and orthonormalize the 128+256 direct
sum there.  Measure the complete affine generator residuals on the same
selected physical-grid seed.  Only after accepted Wannier90 localization may
the existing center-based core-to-buffer materialization redistribute MLWFs.

**Step 3: Preserve receipts and memory bounds**

Publish composition workspace/fingerprint and keep all OA4 occupied receipts.
Density comparisons must use the smoothly composed pre/post occupied
projectors, with strict interior and separately declared boundary tolerances.
Reject rank loss, incomplete coverage, zero workspace after nonidentity work,
non-real Gamma data, or any over-tolerance affine generator.

**Step 4: Focused verification and reviews**

Run EigenExa, construction, checkpoint, DMN/Wannier90, route, obsolete-route,
and Exp-only RT fixtures on MPI 1/2/4/8.  Perform specification and code-quality
reviews for ordering, normalization, buffer lifetime, peak-memory accounting,
MPI determinism, and normal-route isolation.  Resolve every
Critical/Important finding via RED, run `git diff --check`, then run the
clean-first full-feature committed-parent overlay and commit OA5c.

### Task OA5d: Genuine ideal-Si64 GS, MLWF, and V3 acceptance

**Files:**
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`
- Modify tests or production files only for defects exposed by acceptance

#### Task OA5d.1: Replace fixed-center adaptation with cocycle-aware point-cogroup adaptation

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the cocycle RED**

Build a finite affine fixture with a nontrivial translation cocycle.  Start
from a translation-invariant occupied projector and require averaging over one
representative per point coset to produce the same projector as explicit
full-affine averaging.  Require exact rejection when a cocycle entry or
representative product is corrupted.  Extend the route contract so production
must call translation adaptation before point-cogroup adaptation, must pass
`global_point_representatives`, `global_point_cogroup_product`, and
`global_translation_cocycle`, and may use the fixed-center group only after
full-affine acceptance.

Run:

```sh
python3 tests/dg/run_dg_overlapping_wannier_fragment_symmetry_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL on MPI 1/2/4/8 because the current second average uses only the
12 fixed-center operations and the group-average API cannot validate a
translation cocycle.

**Step 2: Implement the minimal cocycle-aware average**

Add a dedicated EigenExa entry point taking translation maps, coset
representative maps, point product, and translation cocycle.  Validate every
representative product against the corresponding representative-plus-
translation action.  Construct the rank-128 averaged projector using the 48
representatives only.  Fill its distributed real Gram matrix from relative
coset/cocycle overlap blocks without retaining transformed real-space orbits.
Reject rank loss, split eigenvalue clusters, non-real Gamma data, zero measured
workspace, overflow, and nonfinite inputs.

In `main_dft.f90`, retain the first translation adaptation, replace the second
fixed-center average by the new point-cogroup average, and publish
`point_cogroup_adapted_occupied`.  Keep the fixed-center group for DMN and
Wannier90 representation generation only.  Do not relax
`dg_ow_symmetry_tolerance` or the final affine-generator gate.

**Step 3: Focused verification**

Run the fragment-symmetry, construction, EigenExa, checkpoint, Wannier90,
route, obsolete-route, and Exp-only RT fixtures on MPI 1/2/4/8.  Require
rank-independent spectra/fingerprints, nonzero bounded workspace, exact
cocycle adverse-case rejection, Gamma-real output, and `git diff --check`.

**Step 4: Specification and code-quality reviews**

Review group multiplication order, pullback action order, cocycle orientation,
normalization by 48 rather than 1536, complete-cluster selection, integer and
MPI-count overflow, error agreement across ranks, memory scaling, and normal
DC isolation.  Resolve every Critical/Important finding with a new RED.

**Step 5: Clean overlay and commit**

Commit the implementation, checkout that commit in the full-feature overlay,
and run:

```sh
cmake --build /tmp/salmon-oa5b-overlay/build --clean-first -j1
```

Rerun focused MPI 1/2/4/8 verification with the clean binary and commit any
review correction separately.

**Step 1: Run immutable genuine Si64**

Use `tests/dg/data/si64_overlapping_wannier_rt/input_gs.in`, its tracked
ideal 64-Si coordinates, and
`samples/exercise_04_bulkSi_gs/Si_rps.dat`.  Run 8 MPI ranks and one OpenMP
thread from a fresh directory with the clean full-feature overlay.  Reject any
run using the older carbon pseudopotential/cell fixture.

**Step 2: Require complete production acceptance**

Require normal DC-SCF/LCFO+EigenExa, occupied rank 128, smoothly composed
complete `s+p` rank 256, final rank 384, affine orders 1536/32/48,
inversion-containing fixed-center group, strict full-affine generator closure,
bounded nonzero workspace, converged symmetry-adapted Wannier90, center-orbit
closure, accepted V3, and bitwise-unchanged checkpoint restart reuse.

**Step 3: Review, overlay, and commit**

Record exact commands, input/pseudopotential/binary/checkpoint hashes,
residuals, density receipts, memory peaks, Wannier90 iterations/spreads, and
restart evidence.  Perform specification and code-quality reviews, resolve
every Critical/Important issue through a fresh RED task, rerun focused MPI
1/2/4/8 and the clean committed-parent overlay, then commit acceptance results.

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
