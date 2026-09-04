# Core-Centered Construction WF Selection Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Select core-centered construction WFs without discarding raw cache data, then admit only a verified selected-WF+PW space to Task 8 local SCF.

**Architecture:** Keep immutable raw construction and selected production catalogs separate. Build the PW complement against selected WFs, recompute the DC seed projection in the actual core metric, and adapt fixed-reference preconditioning to the rectangular selected frame. Resume the existing Task 8 main integration only after the complete small-operator handoff passes.

**Tech Stack:** Fortran 2008, MPI, Wannier90 adapter, BLAS/LAPACK, Python fixture runners, CMake.

---

Read `2026-09-04-hybrid-core-center-selection-design.md` completely first.
This is an amendment to Task 8, not a new task/session/worktree request.
Execute in the user-fixed existing worktree. Preserve all dirty changes and
logs; stage only task-owned hunks. Keep one rank per fragment. Do not rerun
conventional DC or heavy Si64 tests. Use `@test-driven-development`,
`@systematic-debugging` for failures, `@requesting-code-review` at checkpoints,
and `@verification-before-completion` before success claims.

The first execution checkpoint is **Task C1 only**. Later tasks depend on its
center provenance; do not claim selection or Task 8 is implemented after C1.
The preconditioner-frame detail in C4 requires independent design review before
its implementation; do not silently pass a rectangular map to the old API.

### Task C1: Preserve final localized centers in the raw cache

**Checkpoint, 2026-09-04:** Implemented and verified on 2/4/8 MPI ranks;
the release build passes. The raw cache now stores final transform-ordered
`centers_fractional` in `[0,1)`, including canonicalization of a modulo result
rounded to exactly one for a tiny negative coordinate. Cache publication moves
the new array, reuse validates allocation/shape/finiteness/range, and version 2
of the replicated integrity hash includes centers without changing the basis
or transform fingerprint definitions. Reversed/unwrapped stub centers test
alignment with the final transform and values; one-rank center mutation,
missing allocation, wrong extent and NaN are rejected without another W90 call.
The tests caught both a missing publication move and the negative-roundoff
endpoint issue before correction. The unselected core-metric diagnostic still
rejects zero reference norms; C1 does not select WFs or fix that handoff.
C2 and the parent Task 8 production switch remain pending. Existing dirty
code and material validation logs were preserved; no DC run was repeated.

**Files:**
- Modify: `src/gs/dc/dg_hybrid_fragment_wannier.f90`
- Modify: `tests/dg/test_dg_hybrid_fragment_wannier_mpi.f90`
- Test runner: `tests/dg/run_dg_hybrid_fragment_wannier_mpi.py`

**Step 1 — RED test:** Extend the W90 stub with known unequal centers whose
sorting changes the returned column order. Check that cache centers match the
final transform and stored values, not the unsorted run output. Check wrapped
centers, finite 3-by-retained_rank extent, cache reuse without another W90 call,
and collective rejection after one-rank center corruption. Raw retained counts
and every seed reconstruction must be unchanged.

**Step 2 — Run RED:**
`python3 tests/dg/run_dg_hybrid_fragment_wannier_mpi.py`
Expected: missing center cache field or failed center-integrity assertion.

**Step 3 — Minimal implementation:** Add
`real(real64),allocatable :: centers_fractional(:,:)` to the raw cache.
Store the final centers after `apply_dg_w90_gamma_transform`, canonicalized by
the documented fragment-periodic convention. Extend cache allocation/extent,
finiteness, snapshot equality and replicated integrity hashing. Do not change
the basis fingerprint's meaning without explicit versioning; include center
integrity separately if needed. Publish only with the complete valid cache.

**Step 4 — GREEN:** Run the same fixture (2/4/8 rank legacy kernel layouts)
and `cmake --build build-hybrid-release -j4`. Check default fixtures, corruption
failure message and preserved construction call counts. The opt-in unselected
core-metric audit is still expected RED, not fixed by metadata storage.

**Step 5 — Checkpoint:** Review and commit only these task hunks with
`feat(dg): retain final construction Wannier centers`; report C1 complete and
wait for feedback under executing-plans.

### Task C2: Select stable raw columns with periodic core ownership

**Checkpoint, 2026-09-04:** Implemented the geometry classifier and the
single-owner raw-cache selection adapter in `dg_hybrid_fragment_selection`.
Selection returns stable raw column IDs and provenance only; raw WF/buffer
values remain unchanged. The adapter reuses the authoritative cache validator
through coordinate export on the singleton fragment communicator (temporary
Q/seed maps are discarded), validates the actual DC index mapping, and checks
its raw-core-first origin/lattice against the shared grid-aligned core boxes.
Foreign-centered columns are not transferred to another fragment.

The convention uses a common snapped grid coordinate with a 64-epsilon,
geometry-scaled tolerance and periodic half-open boxes. Excessive coordinate
uncertainty is rejected before coordinate division/conversion. The geometry
fingerprint includes the convention; the final receipt also binds the raw
cache, distributed WF integrity and exact rank--fragment inventory.

The new 1/2/4/8-rank fixture passes, including reversed rank ownership,
unequal core widths, translated origins, periodic faces/edges/corners, 3-D
partitions, roundoff boundary probes, empty selections, coincident centers,
phase/permutation transport, corrupt centers/cache/mapping, partition overlap
and rank-disagreeing geometry. Review identified an extreme finite-geometry
overflow path; its trap-enabled RED reproduction now returns a collective
error after correction. The original WF fixture passes on 2/4/8 ranks and the
release build passes. The existing test W90 stubs are reused via a program-only
preprocessor guard; they are not evidence of material localization accuracy.
Review has no remaining Critical/Important issue. C3 seed/PW admission and the
Task 8 main route remain unimplemented; no material/DC calculation was run.
The unrelated CMake fingerprint entry and pre-existing main edits are retained.

**Files:**
- Create: `src/gs/dc/dg_hybrid_fragment_selection.f90`
- Modify: `src/gs/dc/CMakeLists.txt` (only the new module line)
- Create: `tests/dg/test_dg_hybrid_fragment_selection_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_fragment_selection_mpi.py`
- Read: `src/gs/dc/dg_hybrid_fragment_wannier.f90` DC physical mapping

**Step 1 — RED:** Write fixtures for centers inside, outside with nonzero core
tails, on faces/edges/corners, near roundoff boundaries, and across total-cell
and fragment-cell periodic wraps. Include translated origins, unequal core
sizes, reordered rank--fragment mapping, equal-center multiple WFs, zero
selected WFs, NaN centers and inconsistent geometry. Verify selected IDs retain
raw order and input cache values remain bitwise unchanged. Production fixtures
use 1/2/4/8 fragments on the same number of ranks, not split columns.

**Step 2 — Run RED:**
`python3 tests/dg/run_dg_hybrid_fragment_selection_mpi.py`
Expected: missing selection module/API.

**Step 3 — Minimal implementation:** Define a selection receipt containing
raw/selected counts, fragment/generation, raw-cache integrity, geometry,
center convention, selected raw IDs and a nonzero selection fingerprint.
Use a collective API of this shape (exact type definitions belong in this
module, not in main):

```fortran
call select_dg_hybrid_core_wannier(comm_total,fragment_id,cache,&
  fragment_lattice,fragment_raw_origin,total_lattice,total_origin,&
  core_lower,core_extent,raw_grid_shape,core_grid_shape,total_grid_shape,&
  dc_to_total_indices,selection,ok,message)
```

Here `dc_to_total_indices` carries the actual `dc%jxyz_tot` coordinate-index
map, not an inferred origin. Validate aligned geometry against the same mapped
raw core grid before use.
For each center, map to the canonical periodic cell, snap only roundoff-close
shared boundaries and apply the half-open ownership rule. Select iff its owner
equals its construction fragment. Return IDs and receipts only: no W90 call,
foreign-WF redistribution, raw-cache mutation or nonzero-tail-based exception.
No partial selection output after any collective failure.

**Step 4 — GREEN:** Run the new fixture and the original fragment-Wannier
fixture; build. Require identical physical selected spans under raw column
phase/permutation with centers transported consistently. General rotations
mixing selected/excluded WFs are deliberately not a selector invariance test.

**Step 5 — Commit:** Only C2 files/hunks, message
`feat(dg): select core-centered construction Wannier columns`.

### Task C3: Selected PW complement, core metric and DC projection

**Partial checkpoint, 2026-09-04 (selected-value export only):** Added
`export_dg_hybrid_selected_wannier`. It revalidates the immutable raw cache,
checks selection metadata/ordered raw IDs and the current exact rank--fragment
inventory, then exports selected columns without removing any buffer samples.
Selection receipt hash version 2 additionally seals generation, counts,
convention and center owners. Empty selection is still valid in the selector,
but this downstream export collectively rejects unsupported PW-only fragments.
No raw seed coefficient slicing, occupation change or W90 rerun is introduced.

The new API initially failed to link (expected missing implementation); the
1/2/4/8-rank fixture now passes full-value/order checks, single-rank receipt or
cache corruption rejection, missing IDs, empty selection and raw immutability.
The existing 2/4/8-rank raw-Wannier fixture and release build pass. Independent
review found no Critical/Important issue; temporary validation maps are freed
before copying values. This is only the first C3 connection component: active
global IDs, selected-only PW projection, actual core Gram solve, reconstruction
diagnostics and the fail-only post-initializer density gate remain pending.
The production main route remains unchanged; no DC/material run was repeated.

**Files:**
- Modify: `src/gs/dc/dg_hybrid_fragment_selection.f90`
- Modify: `src/gs/dc/dg_hybrid_projected_fragment_pipeline.f90`
- Modify: `tests/dg/test_dg_hybrid_fragment_selection_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_fragment_selection_mpi.py`
- Modify: `tests/dg/test_dg_hybrid_projected_fragment_pipeline_mpi.f90`
- Reuse: `src/gs/dc/dg_hybrid_broken_volume.f90`

**Step 1 — RED:** Use a nontrivial localized transform with an excluded
external-centered WF having finite core tail. Demonstrate that slicing the
old seed coefficients fails, then require projection onto selected WF+PW to
match an independent physical-space least-squares oracle. Include a case
where projection against the raw rather than selected WF union incorrectly
eliminates a needed PW component. Include insufficient user cutoff, residual
metric dependence, invalid selection fingerprint and selected-empty cases.
Separately include an exactly representable DC seed with core norm 1/2: the
current initializer changes its density if occupations are left unchanged.
Require the handoff to reject this specific post-initializer density mismatch,
not to report a PW span failure or normalize the density silently.

**Step 2 — Run RED:** Run the selection and projected-fragment MPI runners.
Expected: missing selected-catalog handoff or failed seed/span assertion.

**Step 3 — Minimal implementation:** Form selected WF views from receipt IDs
without dropping any buffer rows. Build global active-column IDs distinct from
stable raw IDs and pass only selected WFs to the generalized PW projection.
Retain the union-complete map separately. Reject unsupported PW-only input
explicitly rather than silently restoring excluded WFs.

Build the actual core Gram S and certify its resolved positive rank; use the
existing tolerance policy, not diagonal positivity alone. With
`rhs=B_core^dagger W_core Psi_core`, solve `S C=rhs`. Preserve raw seed spectra
and occupations as reference metadata; do not silently infer physical state
counts from selected basis count. Record orbital, occupied density/electron
and required support reconstruction defects independently. Emit a collective
insufficient-span diagnostic with selected count and cutoff on failure.

**Step 4 — GREEN:** Run both focused runners, the raw-Wannier runner and build.
Certify projection density before initialization, then compare after the
existing initializer with the actual occupations. This increment is fail-only
admission: if post-initializer density/electron tolerance fails, reject without
publishing a solver state; no occupation transformation is implemented here.
Require all-retained selection with an admissible core metric/seed to recover
the previous physical result and no extra W90 invocation. If representative
small DC inputs require general density-preserving occupation transformations,
pause at this checkpoint for explicit design before C5/C6. Positive toy seeds
alone must not be advertised as a general production handoff.

**Step 5 — Commit:** C3 hunks only, message
`feat(dg): project DC seeds into selected WF plus PW space`.

### Task C4: Fixed projected-reference preconditioner

**Files:**
- Modify: `src/gs/dc/dg_hybrid_fragment_preconditioner.f90`
- Modify: `tests/dg/test_dg_hybrid_fragment_preconditioner_mpi.f90`
- Modify: `src/gs/dc/dg_hybrid_fragment_selection.f90`
- Test: `tests/dg/run_dg_hybrid_fragment_preconditioner_mpi.py`

**Step 1 — Review design:** Independently review the rectangular frame formula
in the design amendment. Resolve any change in safety or invariance semantics
before writing code. A requirement for a different preconditioner is a design
checkpoint, not authority to substitute identity.

**Step 2 — RED:** Define explicit rectangular Q fixtures with `Q Q^dagger=I`,
nonidentity core S, nontrivial shifts, and a physical-space expected action.
Test selected-space full unitary covariance, exact-zero reference columns,
nonzero columns with unresolved S norm, rank-deficient frame, changed selection
fingerprint, stale epochs and the square-limit equivalence to the existing API.
Run the preconditioner runner and observe the new frame-entry failure.

**Step 3 — Minimal implementation:** Add an explicit frame-preparation API;
do not weaken the existing square-unitary entry. Reference column count may
exceed active row count. Omit only roundoff-zero coefficient columns, recheck
the identity resolution and keep all finite/metric/shift safeguards. Form
`h_a=q_a^dagger H q_a`, `s_a=q_a^dagger S q_a` and apply the existing signed
regularization in `sum_a q_a D_a^-1 q_a^dagger r`. Bind selection provenance.
Never rerun localization or change the frame during the density loop.

**Step 4 — GREEN:** Run preconditioner, selection, raw-Wannier and bounded
subspace fixtures plus build. Verify unchanged square-action results and
physical, not just coefficient-space, covariance after selection.

**Step 5 — Commit:** C4 hunks only, message
`feat(dg): support fixed projected frames after WF selection`.

### Task C5: Actual small DG handoff and the known counterexample

**Files:**
- Modify: `tests/dg/test_dg_hybrid_fragment_wannier_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_fragment_wannier_mpi.py`
- Modify: `tests/dg/test_dg_hybrid_fragment_selection_mpi.f90`
- Reuse production volume, interior/face, nonlocal and payload modules.

**Step 1 — RED:** Keep the unselected buffer-only counterexample, but assert
the specific zero-reference-norm rejection rather than any nonzero exit. Add
a center-selected positive case with physically consistent stub centers and
localized values. A separate tail case must fail below a sufficient PW cutoff
and pass at a user-explicit sufficient cutoff. Use actual volume, SIPG and
nonlocal assembly for acceptance, with independent small quadrature oracles.
Zero gradients/toy H alone cannot certify this stage.

**Step 2 — Run RED:** Run both fixtures; verify the selected path, not a parser
or missing unrelated library, is responsible for the expected failure.

**Step 3 — Connect:** Pass selected receipt, core S, projected seeds and frame
to the real payload/self-block/bounded-CG APIs. Use the selected basis for all
volume, derivative, face, nonlocal and density calculations. Require safe
metric orthogonality, independent physical density/core norms, no loss of
required boundary/projector data, and a shared maximum-three-update budget
across extensions in one density epoch. No new numerical bypass is allowed.

**Step 4 — GREEN:** Run the full relevant small MPI set on the production
one-rank-per-fragment topology and build. Confirm the negative rejection reason
and positive selected-path physics separately. Raw cache W90 counts stay one.

**Step 5 — Checkpoint:** Review C1--C5 together; commit only current hunks,
then report readiness to resume Task 8, not completion of Task 8.

### Task C6: Resume the separated Task 8 production route

**Files:**
- Modify: `src/gs/main_dft.f90` (task hunks only; preserve prior dirty code)
- Modify: `tests/dg/check_dg_hybrid_fragment_wannier_route.py`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Update: `docs/plans/2026-09-02-hybrid-divided-one-shot-lcfo.md`

**Step 1 — RED route contract:** Require direct raw DC construction, verified
center selection before selected-WF PW projection, core-metric/seed admission,
frame preparation and bounded updates. Forbid using raw retained_rank as the
production active dimension or padding raw coefficients without projection.
Keep all existing no-preliminary-LCFO/global-W90/dense-callback guards.

**Step 2 — Run RED:**
`python3 tests/dg/check_dg_hybrid_fragment_wannier_route.py`.
It remains intentionally RED until the actual separate main path is connected.

**Step 3 — Implementation:** Resume Steps 3--4 of parent Task 8 with this
amended catalog/seed/frame flow. Change all shape/ownership/fingerprint users
together; do not promote a partially connected entry based only on call names.
Keep construction/selection outside SCF and LCFO after converged density.

**Step 4 — GREEN:** Run parent Task 8 focused tests and build. Main completion
also requires the operator integration results from C5 and the original
occupation/density/budget gates, not just source-token checks. Material and RT
validation remain the later parent tasks using existing DC files.

**Step 5 — Checkpoint:** Update parent progress accurately, review and commit
only this task's hunks; report remaining parent tasks explicitly.
