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

**Second partial checkpoint, 2026-09-04 (selected PW connection and core solve):**
`prepare_dg_hybrid_selected_catalog` now exports validated selected values with
compact fragment-major active IDs, a separate replicated raw-column provenance
directory and a shared fingerprint binding every selection receipt and the
exact rank--fragment inventory. Unequal local counts are supported; generation
disagreement and invalid/empty selection cannot publish a catalog. This catalog
feeds the existing local projected-fragment pipeline without adding a dependency
from that generic pipeline back to Wannier construction.

The 2/4/8-rank integration fixture uses the real raw-cache constructor (only W90
is stubbed), center selector, selected export and PW projection. An excluded WF
has a finite core tail: projection against the raw union removes a needed PW
component, whereas the selected union retains it. All selected buffer values
and raw/active ID mappings remain intact. The centers/values are a controlled
algebraic fixture, not evidence of material localization accuracy.

`project_dg_hybrid_core_seeds` builds the weighted core normal equations and
uses the existing Hermitian metric rank policy. It rejects any rank loss before
publishing coefficients, rather than adopting a compressed local basis. Its
local report separates maximum relative orbital residual, relative weighted
density L1 defect and absolute electron-number defect. Insufficient span reports
selected count, supplied cutoff and measured defects; no cutoff or occupation
is modified. The cutoff is diagnostic input, not a new PW catalog generator.

Tests compare a nonorthogonal coordinate example to an independent coefficient
oracle and to the actual broken-volume unit-potential metric. Raw DC orbitals
are reconstructed from the full immutable cache map before projection; their
spectra and occupations are preserved. Unequal per-rank basis/seed dimensions,
missing PW span, dependent metric, finite overflow and differing controls are
covered. The selection runner passes on 1/2/4/8 ranks (PW/core integration on
2/4/8), both original 2/4/8-rank WF/pipeline runners pass, and release builds.

C3 remains incomplete: production admission still needs the bound raw-DC seed
adapter, required boundary/derivative/projector support diagnostics, and the
actual initializer's fail-only density/electron gate (including core norm 1/2).
The probe above does not authorize publishing an accepted solver state or
claim a general density-preserving DC handoff. C4/C5/C6 and main are unchanged;
no material/DC calculation was rerun and all pre-existing dirty data is retained.

**Third partial checkpoint, 2026-09-04 (actual initializer density gate):** Added
`initialize_dg_hybrid_fragment_density_checked` in the subspace module, preserving
the legacy initializer API. The explicit single-owner adapter validates a
bijective rank--fragment map, initializes a temporary state on MPI_COMM_SELF,
and computes physical core density with the returned selected seed IDs and
their original occupations. The caller's reference density is never rescaled.
Relative weighted density L1 and absolute electron defects are returned
separately. Any rank's failure preserves every caller state and publishes no
selected-seed output; state vectors are moved only after collective admission.
The callback must be fragment-local, without total-communicator collectives.
Tolerances, guard count and optional cutoff presence/value are shared controls;
local seed spectra, occupations and dimensions may differ.

The integrated 2/4/8-rank test first passes an exact core projection with seed
norm 1/2, then verifies that the real initializer doubles density and is rejected
specifically as a post-initializer density mismatch, not insufficient PW span.
A unit-core-norm positive case passes. Fractional occupation/energy reordering,
equal-electron-number density redistribution, one-rank failure rollback, metric
callback failure, NaN reference and differing controls are also covered. Review
prompted the guard/cutoff agreement tests; the guard mismatch was reproduced RED
before adding those checks. The selection runner passes on 1/2/4/8 ranks, the
legacy subspace runner passes on 1/2/4/8, the raw-Wannier runner passes on 2/4/8,
and release builds. Independent review has no remaining Critical/Important issue.

This is fail-only admission, not an occupation transformation or a completed
production DC handoff. C3 still needs authoritative raw-reference/selected-basis
binding and required boundary/derivative/projector reconstruction diagnostics.
The half-norm example is controlled test data, not a representative material
calculation triggering the separate occupation-design decision. No C4/C5/C6 or
main switch was attempted; old dirty changes and verification logs are retained.

**Fourth partial checkpoint, 2026-09-04 (bound raw DC reference export):**
The selection receipt now retains the actual mapper's periodic physical grid
IDs and core row slots; version 3 of its integrity hash includes their lengths
and contents. Selected exports reject absent, invalid or modified row mappings
before using them. This preserves arbitrary raw storage order rather than
assuming that core samples occupy a contiguous prefix of the stored array.

`export_dg_hybrid_dc_reference` validates both cache and selection, reconstructs
buffer orbitals from **all** raw WFs and the original DC coefficient map, and
extracts core orbitals using the sealed row slots. It returns the unchanged DC
energies/occupations, physical row IDs, generation and local selection binding.
No selected-column slicing, localization rerun or density normalization occurs.
All outputs remain unpublished on any rank's validation failure.

Tests cover unequal core sizes, periodic mapping, deliberately interleaved
core/buffer storage, corrupt row slots/IDs/seed coefficients/occupations, and
immutable raw reconstruction. The small integration test now uses this exported
reference for the core solve and density-checked initializer. The selection
runner passes on 1/2/4/8 ranks, raw-Wannier and projected-pipeline regressions
pass on 2/4/8, and release builds. Independent review found no Critical/Important
issue in the export binding. This is not a persisted reference-file format.

C3 remains open for required boundary/derivative/projector support diagnostics
and the final combined admission path that verifies its mutable selected basis
and reference together. The generic numerical kernels do not independently
validate an externally modified reference object. Production main and C4--C6
remain unchanged; no conventional DC/material run or worktree creation occurred.

**Fifth partial checkpoint, 2026-09-04 (support reconstruction numerical gate):**
Added `s_dg_hybrid_support_samples` and `check_dg_hybrid_seed_support` to the
projected pipeline module. Boundary samples, derivative samples and nonlocal
projector overlaps are checked separately for every seed using the already
computed coefficients. Each reported defect is the maximum absolute weighted
L2 orbital error in that channel's physical units, following the raw seed span
norm convention. Three explicit caller tolerances and cutoff must agree across
ranks; there is no occupation-based waiver or implicit tolerance relaxation.

The numerical API takes independently supplied required sample counts and
rejects missing arrays, mismatched extents, nonpositive/nonfinite weights,
nonfinite values or arithmetic overflow. Allocated zero-row evidence is allowed
only for a declared zero required count; counts may vary by fragment. All three
defects remain available on a measured tolerance failure. No coefficient,
occupation, cutoff or solver state is modified or published by this check.

The integration fixture evaluates explicit boundary sampling, small difference
stencils and a normalized projector functional on real cached raw DC orbitals
and the selected WF+PW values. Exact DC reconstruction passes. A core-exact
excluded-WF probe fails support admission with independently predicted errors
sqrt(1/2), sqrt(1/2)/2 and sqrt(1/6); its buffer tail cannot be waived by the core
density result. Single-channel one-rank defects, quadrature scaling, missing
samples, NaN, finite overflow, differing controls and a declared empty local
projector inventory are covered. Selection tests pass on 1/2/4/8 ranks (support
integration on 2/4/8), existing raw-Wannier/pipeline regressions pass on 2/4/8,
and release builds. Independent review found no Critical/Important issue.

This is numerical admission of supplied evidence, not certification that the
production operator inventory is complete. The final adapter must still bind
actual sample IDs/order and required counts to the selected basis, raw reference
and operator provenance before combined state admission. Production SIPG and
nonlocal operator validation remains C5; these explicit stencils are not a
substitute. C3 remains open, main/C4--C6 are unchanged, and no DC/material run
was repeated. All pre-existing dirty files and verification logs were retained.

**Sixth partial checkpoint, 2026-09-04 (combined numerical admission):**
Added `dg_hybrid_fragment_admission` as a separate adapter. It rebuilds the
selected catalog and raw DC reference, verifies the projected basis payload
receipt, and binds the ordered physical core IDs and quadrature weights to
those used by the PW producer. Changed weights are rejected as an input
binding mismatch, not misreported as a physical density failure.

Frozen sparse support manifests bind channel, fragment, generation, required
sample identities/order, physical grid IDs, coefficients and weights. The
adapter evaluates each functional on both the selected basis and raw DC
reference; callers cannot substitute precomputed support values or density.
Core projection, support checks and density-checked initialization are combined
without publishing an intermediate solver state. Independent review identified
that initialization can change support even when its density error is allowed.
The initializer now writes a temporary state; all three support channels are
checked again against the selected raw seeds before final publication.

The regression uses a separately constructed, extended-domain-normalized seed
with core norm 1/2. Its exact core and support projection passes, but subsequent
normalization fails the density gate. With only the density tolerance relaxed,
the old combined entry incorrectly accepted it (RED); the final support gate
now rejects it (GREEN). Changed basis payloads, stale operator fingerprints,
misordered manifests and changed quadrature also reject without state changes.
No additional W90 calls occur during admission of either cached generation.

Selection tests pass on 1/2/4/8 ranks; projected pipeline and raw-Wannier tests
pass on 2/4/8; fragment-subspace tests pass on 1/2/4/8. Release build succeeds.
Independent re-review has no remaining Critical/Important issue. This remains
a numerical checkpoint: actual production operator inventory adapters and
SIPG/nonlocal acceptance are C5, and the main route remains C6. C3 is not being
declared a general physical DC handoff. No DC/material calculation was repeated;
all unrelated dirty changes and existing verification logs were preserved.

**Files:**
- Modify: `src/gs/dc/dg_hybrid_fragment_selection.f90`
- Modify: `src/gs/dc/dg_hybrid_projected_fragment_pipeline.f90`
- Modify: `src/gs/dc/dg_hybrid_fragment_subspace.f90` (density-checked initializer adapter)
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

**Approved continuation — first increment:** Retain signed denominators and
add fail-only cancellation detection, as approved by the user after the review
below. First implement and test the numerical gate separately; wire it into
the new rectangular entry in the next increment, without changing square API
semantics. For each state, with signed-scaled reference amplitudes y, compare
`||Q y||_2` with `tau * sum_a ||q_a||_2 |y_a|`, where
`tau=max(tolerance,64*epsilon_machine*max(active_count,reference_count))`.
Use collective Euclidean norms (invariant under a unitary active-space change),
scale amplitudes before squaring, and accept an identically zero sum of term
norms. A resolved term sum but unresolved result is an explicit collective
failure. This diagnoses cancellation for the applied residual; it does not
certify nonsingularity for every possible residual.

RED/GREEN for this first increment: extend
`tests/dg/test_dg_hybrid_fragment_preconditioner_mpi.f90` with the reviewed
three-column example, a noncancelling residual, zero amplitudes, unitary
transport, stale/differing controls and finite/nonfinite inputs. Add a
numerical gate in `src/gs/dc/dg_hybrid_fragment_preconditioner.f90`, run
`python3 tests/dg/run_dg_hybrid_fragment_preconditioner_mpi.py` on 1/2/4/8 ranks,
then review and checkpoint. The rectangular producer, provenance binding,
zero-column omission rule and floor-dimension policy remain subsequent C4
work, not implied complete by this gate.

**Numerical gate checkpoint, 2026-09-04:** The missing gate produced the
expected RED at link time; the implemented gate passes the full preconditioner
runner on 1/2/4/8 ranks, including distributed nonzero row contributions,
complex unitary transport, exact and near cancellation, resolved action,
zero amplitudes, amplitude scales 1e200 and 1e-200, differing controls and
replicated amplitudes, and nonfinite input. Release build succeeds. Independent
review found no Critical/Important issue; its distributed-sum test suggestion
was added and passed. This is a stand-alone numerical check only: no existing
square action or production call path was modified. The next C4 increment
must build and certify the rectangular frame and invoke this gate before
publishing its applied result.

**Design checkpoint, 2026-09-04 — implementation paused:** Independent review
found that the proposed rectangular signed-frame action can annihilate a
nonzero residual despite `Q Q^dagger=I`, positive S and resolved denominators.
For `Q=sqrt(2/3)*[[1,-1/2,-1/2],[0,sqrt(3)/2,-sqrt(3)/2]]`,
`H=diag(1,-1)`, `S=I`, and shift zero, the reference H diagonals are
`[2/3,-1/3,-1/3]` and the resulting action is `diag(0,-3)`.
An independent arithmetic check gives the first diagonal as 4.44e-16;
the signed denominator floor does not address this cancellation. Extending
Q by a scalar identity, taking `H=[[1,0,1],[0,-1,0],[1,0,0]]`, S=I and x=e3
gives an actual Rayleigh residual e1 at shift zero, also annihilated.

Left-unitary covariance and the square limit hold, but do not imply a usable
rectangular preconditioner. Before coding, obtain a decision between explicit
fail-only cancellation detection and a revised action/sign policy. No identity
fallback or unsigned-denominator replacement is authorized by this review.
Also specify the roundoff-zero column threshold and the active/reference
dimension used in the denominator floor. Existing square preconditioner tests
still pass on 1/2/4/8 ranks. No implementation or production-route change was
made at this checkpoint, and no DC/material run was repeated.

**Second C4 partial checkpoint, 2026-09-04 — rectangular numerical entry:**
Added `prepare_dg_hybrid_frame_preconditioner` for a single-owner fragment
communicator (production MPI_COMM_SELF). Multi-rank fragment communicators
are rejected rather than enabling column distribution. The supplied row IDs,
Hermitian H/S and finite inputs are validated; only coefficient columns with
Euclidean norm <= `64*epsilon_machine*max(n,m_input)` are omitted. The retained
frame must still satisfy `Q Q^dagger=I` within the existing tolerance. Nonzero
columns with unresolved reference metric norms are rejected, not removed.
The signed denominator roundoff uses `max(n,m_retained)`; the old square entry
therefore retains its previous dimension and sign convention.

The private cache binds the original Q payload (including omitted columns),
selection fingerprint, operator key and layout. Applying a rectangular cache
requires the matching selection fingerprint. After the signed action is
formed, the cancellation gate runs before output publication. Failed rebuilds
preserve the prior cache. This routine checks reference metric norms, not full
S rank: the separate C3 core-metric admission remains a production prerequisite.

Review found that reusing identity-reference validation also introduced a
spurious coordinate-dependent metric-norm check. A regression with
`S=diag(1e-9,1)` and Hadamard Q demonstrated RED versus the old square API.
Matrix validation is now separated privately from actual-reference metric
validation; no public skip-validation option exists. The regression is GREEN.

Fresh preconditioner tests pass on 1/2/4/8 ranks, including a nonidentity-metric
oracle, complex unitary covariance, square-limit equivalence, zero columns,
unresolved nonzero columns, rank-deficient frames, stale/missing selection,
operator epochs, optional-argument disagreement, rollback and cancellation
without publication. Selection/subspace tests pass on 1/2/4/8, raw-Wannier on
2/4/8, and release build succeeds. Independent re-review has no remaining
Critical/Important issue. The authoritative raw-U/selection-to-frame exporter
and its integration with selected admission remain next; C4 and Task 8 are
not complete. Existing dirty files/logs were retained and no DC run repeated.

**Third C4 checkpoint, 2026-09-04 — cached reference connection:** Added
`export_dg_hybrid_selected_frame` to revalidate the raw cache and selection,
then return `E^dagger U^dagger` using stable raw column IDs, not a prefix or a
new localization. Its reference fingerprint binds selection, transform and
exported coordinates. `export_dg_hybrid_selected_basis_frame` rebuilds the
selected catalog, validates the actual projected-basis receipt, and verifies
WF values, active IDs, sectors, generation and physical row layout before
appending the independent PW identity block. The augmented reference receipt
also binds the complete selected-WF/PW payload; failures publish no frame.

The fixture connects the real cached transform, selected catalog and PW
producer to the real single-owner rectangular preconditioner. A physical-grid
sum over projected reference functions agrees with its action. Unequal local
counts and non-prefix selected raw IDs are tested separately. Corrupt selection
metadata, changed raw transforms and stale projected receipts fail collectively;
W90 setup/run counts do not increase. The initial missing-exporter RED is now
GREEN. Selection/preconditioner/subspace tests pass on 1/2/4/8 ranks, raw-Wannier
on 2/4/8, and release build succeeds. Independent review has no remaining
Critical/Important issue.

C4's numerical reference construction and application are now connected. This
does not certify an actual DG operator handoff: the fixture's explicit local
potential is an arithmetic oracle, not a substitute for C5 SIPG/nonlocal tests.
Next execute C5 with the separate C3 core/support/state admission gates and
the bounded local update. Main remains unchanged until C6; Task 8 is not
complete. All pre-existing dirty changes and verification logs were retained,
and no material/DC calculation was repeated.

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

**Design pause, 2026-09-04 — saved representative DC evidence:** A read-only
audit of the existing eight-rank Si64 conventional DC checkpoint now triggers
the explicit C3 occupation-design stop condition. Occupied orbitals have unit
extended-grid norm but core norms 0.0147253..0.4806794. Original occupations
give about 32 core electrons per fragment; retaining those occupations after
core orthonormalization would instead give about 168.0855. The core overlap
also has sizable off-diagonal entries. This is not a PW-cutoff deficiency and
cannot be resolved by the existing fail-only admission policy.

See `../notes/2026-09-04-dc-core-density-handoff-audit.md` for checkpoint identity,
binary layout, physical core mapping, manifest-weight cross-check and results.
No DC/W90/SCF run was repeated. This is an independently computable necessary
condition, not a claim that the main route was connected and run. Pause C5/C6
promotion for user approval of density-preserving occupation transport or a
different explicit starting-density policy. All prior small numerical tests
remain useful but cannot establish acceptance of this representative DC state.

**First assembly checkpoint, 2026-09-04:** The old unselected audit was first
rerun and reproduced its RED at the preconditioner, with independently computed
zero reference norms in both fragments. It is now a permanent expected-negative
test: require exactly `nonpositive or unresolved reference metric norm`, a
zero output fingerprint and an independently verified zero norm; return before
initialization. The normal raw-Wannier runner also executes this two-rank audit
and requires a success marker from each fragment. Arbitrary failures cannot
pass this negative case.

The selected-WF/PW fixture now calls actual broken-volume, SIPG and nonlocal
row assemblers, freezes their single-owner payload and extracts the local
self-block. Core derivatives are explicit nonzero one-sided/central differences;
face values extrapolate both neighboring cores to their common midpoint with
a consistent +x derivative. Each periodic interface is assembled once. A
normalized complex two-point projector on each interface spans both fragments.
Independent quadrature/trace/projector sums check the matrices, including
nonzero cross-fragment blocks and retention of every term in the self-block.
This is a test-only small discretization, not certification of the eventual
production trace/stencil/projector support inventory.

Selection tests pass on 1/2/4/8 ranks (the assembly case on 2/4/8, exactly one
rank per fragment). Raw-Wannier tests pass on 2/4/8 plus the expected-negative
audit on 2; existing SIPG tests pass on 1/2/4 and divided-operator tests on
1/2/4/8. Release build succeeds. Independent review has no Critical/Important
issue for this assembly checkpoint. No production source was changed, no DC
calculation was repeated, and existing dirty changes/logs were retained.

C5 remains open: connect independently complete operator support manifests,
exercise an explicit insufficient/sufficient PW cutoff tail case, and pass the
same C3-admitted state through these actual operators and the bounded CG budget.
Do not infer those properties from separate passing assembly and solver tests.
Main/C6 and general material acceptance remain untouched.

**Second integration checkpoint, 2026-09-04:** The assembled volume/SIPG/
nonlocal self-block now advances the same state published by combined C3
admission, without reconstructing the state or rewriting its basis/metric
fingerprints. Frozen support manifests cover both faces, all four derivative
rows and both projectors touching each core, with coefficients matching this
small fixture's assembly maps. A missing neighboring projector contribution
is rejected collectively without publishing a nonlocal matrix.

The authoritative cached selected frame feeds the real rectangular
preconditioner and bounded epoch updater. Two calls share one density epoch;
their cumulative updates stay within three and the returned remaining budget
matches that cumulative count. The initial residual is nonzero and the final
occupied physical subspace changes, excluding a no-op or only an internal
rotation. After each call, metric orthogonality, reconstructed core density,
electron count and physical core norms agree with independent grid sums.
This does not force budget exhaustion in the tiny fixture or certify the
post-exhaustion measurement-only branch; that stronger claim is not made.

Fresh selection/preconditioner/subspace tests pass on 1/2/4/8 ranks (this
integration on 2/4/8), raw-Wannier tests pass on 2/4/8 with the exact core-null
negative on 2, and release build succeeds. Independent review has no new
Critical/Important issue. W90 call-count checks remain unchanged and pass.
Only tests and this progress record changed; no DC/material run was repeated.
All unrelated dirty changes and verification logs were preserved.

C5 remains open for the explicit insufficient/sufficient user-cutoff tail
case and certification of production support-provider completeness. The
matching hand-written fixture manifests are not that production certificate.
The same-state bounded-update connection is now tested, but C5/Task 8 and the
main/C6 production switch are not complete.

**Cutoff projection diagnostic checkpoint, 2026-09-04:** Added a separate
analytic tail probe, not a raw-cache localization or combined DC admission
fixture. It uses eight core points, a retained core delta and an unchanged
extended reference with a sine core tail (core norm 1/2). The real reciprocal
catalog selects G=0 at explicit cutoff zero and G=0,+/-G at the explicit first
shell energy. Those catalogs feed the real selected-WF/PW projection pipeline.
Both core metrics retain full rank. Low cutoff fails with a measured orbital
residual matching an independent arithmetic-mean fit; high cutoff recovers
the core tail and its density/electron count without changing the reference,
occupations or requesting another W90 run. No initializer is called: this
half-core-norm reference must not be advertised as an admitted DC state.

The new trap-enabled test exposed an existing overflow in the PW phase safety
guard itself: `(huge/4)/g_scale` overflows for small nonzero g_scale. It passed
on two ranks but reproduced SIGILL on four ranks, including after moving the
test before the old integration fixture. The minimal production fix evaluates
that quotient only for g_scale>1/4. At or below 1/4 every finite coordinate
product is already bounded by huge/4, so the three-term phase sum remains safe.
No tolerance, floating-point trap or unsafe-large-phase rejection was disabled.

The new regression passes on 2/4/8 ranks; the full selection runner passes on
1/2/4/8. Windowed-PW (including unsafe-large-phase rejection) and reciprocal
catalog tests pass on 1/2/4/8, projected pipeline tests on 2/4/8, and release
build succeeds. Independent review found no Critical/Important issue in the
guard fix or diagnostic scope. Existing dirty changes/logs were preserved.

This completes only the cutoff-to-core-projection diagnostic component.
C5 still requires cutoff comparison with bound raw DC states and actual
operator support admission, plus production support-provider completeness.
The analytic probe does not replace those gates or authorize C6/main changes.

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
