# Windowed-PW/Wannier Hybrid Implementation Plan

> Historical/removed references below identify the superseded monolithic WPW
> implementation only; they are not instructions to restore that route.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a symmetry-complete windowed-PW complement to retained Wanniers, retain a sparse PW metric, and propagate the hybrid basis in the length gauge without any global `N_W x N_P` projection matrix.

**Architecture:** Build the feature as new, narrow primitives around the current overlapping-Wannier and full-cell `hpsi` infrastructure. Reuse numerical ideas from the deleted WPW history only after reducing them to row-owned/window-local interfaces; do not restore its production context, fallback routes, or globally retained projection arrays. Establish correctness with full-cell oracles before enabling sparse neighbor truncation or production RT.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, SALMON full-grid `hpsi`, current overlapping-Wannier checkpoint/RT modules, Python MPI fixture runners.

---

## Working Rules And Review Gates

- Use `@superpowers:test-driven-development` for every implementation task.
- Run production-like tests with `OMP_NUM_THREADS=1`; do not change the requested MPI rank count in production checks.
- Preserve unrelated dirty files. Stage only files named by the current task.
- Every allocation and shape-dependent collective needs checked wide extents, MPI rank agreement, allocation consensus, and failure cleanup.
- After Tasks 3, 6, 9, and 12, stop for two reviews:
  1. specification/numerical review;
  2. quality/MPI/memory/scaling review.
- Resolve all Critical and Important review findings and rerun the focused MPI 1/2/4/8 suites before continuing.
- Use `@superpowers:verification-before-completion` before any completion claim.

### Task 1: Add A Route Guard And Minimal Hybrid Catalog Type

**Files:**
- Create: `src/common/dg_hybrid_windowed_pw_types.f90`
- Modify: `CMakeLists.txt`
- Create: `tests/dg/check_dg_hybrid_windowed_pw_route.py`

**Step 1: Write the failing route test**

Require new `dg_hybrid_windowed_pw_*` modules and forbid production references to deleted monolithic types such as `s_dg_wpw_s_orthogonal_complement`, `a_owned_w_global_p`, and `dg_wpw_production_context`.

**Step 2: Run the test and verify failure**

Run: `python3 tests/dg/check_dg_hybrid_windowed_pw_route.py`

Expected: FAIL because the new catalog/type module and CMake entry do not exist.

**Step 3: Add the minimal types**

Define only metadata initially:

```fortran
type, public :: s_dg_hybrid_pw_packet
  integer :: fragment_id=0, star_id=0, owner_rank=-1
  integer, allocatable :: g_indices(:)
end type

type, public :: s_dg_hybrid_basis_catalog
  logical :: valid=.false.
  integer(8) :: wannier_fingerprint=0_8, window_fingerprint=0_8
  integer(8) :: packet_fingerprint=0_8, catalog_fingerprint=0_8
  integer, allocatable :: accepted_wannier_blocks(:), rejected_wannier_blocks(:)
  type(s_dg_hybrid_pw_packet), allocatable :: packets(:)
end type
```

Do not add operator storage yet.

**Step 4: Run route and diff checks**

Run:

```bash
python3 tests/dg/check_dg_hybrid_windowed_pw_route.py
git diff --check
```

Expected: PASS.

**Step 5: Commit**

```bash
git add CMakeLists.txt src/common/dg_hybrid_windowed_pw_types.f90 tests/dg/check_dg_hybrid_windowed_pw_route.py
git commit -m "feat: add windowed PW hybrid catalog"
```

### Task 2: Select Wannier Symmetry Blocks Atomically

**Files:**
- Create: `src/gs/dc/dg_hybrid_wannier_selection.f90`
- Create: `tests/dg/test_dg_hybrid_wannier_selection_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_wannier_selection_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing MPI fixtures**

Cover:

- accept/reject of a complete translation/point-cogroup block;
- rejection when only one member appears localized;
- conjugate-pair atomicity;
- repeated-center blocks;
- rank-disagreeing thresholds/metadata;
- nonfinite spread/localization receipt;
- permutation and MPI-decomposition invariant selection fingerprint.

**Step 2: Run on MPI 1/2/4/8 and verify failure**

Run: `python3 tests/dg/run_dg_hybrid_wannier_selection_mpi.py`

Expected: FAIL because the selector is absent.

**Step 3: Implement the collective selector**

The selector consumes canonical block membership and localization receipts. It
returns accepted and rejected block IDs, target complement rank per block, and a
fingerprint. It may not inspect or reject isolated columns.

**Step 4: Rerun MPI and route tests**

Expected: PASS on 1/2/4/8 with identical fingerprint.

**Step 5: Commit**

```bash
git add CMakeLists.txt src/gs/dc/dg_hybrid_wannier_selection.f90 tests/dg/test_dg_hybrid_wannier_selection_mpi.f90 tests/dg/run_dg_hybrid_wannier_selection_mpi.py
git commit -m "feat: select localized Wannier symmetry blocks"
```

### Task 3: Build Covariant Partition Windows And Complete G Stars

**Files:**
- Create: `src/common/dg_hybrid_windowed_pw_basis.f90`
- Create: `tests/dg/test_dg_hybrid_windowed_pw_basis_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_windowed_pw_basis_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing fixtures**

Test a periodic Z2xZ2 fragment grid and a nonsymmorphic/permuted fragment catalog.
Verify `sum_f chi_f^2=1`, window covariance, complete G-star membership, conjugate
pairing, deterministic packet order, row-owned materialization, and bounded tile
workspace. Add corrupt window, incomplete star, duplicate row, rank disagreement,
and finite-huge adverse cases.

**Step 2: Run MPI 1/2/4/8 and verify failure**

Run: `python3 tests/dg/run_dg_hybrid_windowed_pw_basis_mpi.py`

**Step 3: Implement only window/packet construction**

Use historical `dg_wpw_windows.f90` and `dg_wpw_g_modes.f90` as mathematical
references via `git show 4c0e7efb^:<path>`, but write a new row-owned API. Generate
PW values in bounded tiles; never retain every packet on every grid row.

**Step 4: Verify**

Run the focused suite, route checker, existing construction MPI suite, and
`git diff --check`.

**Step 5: Review Gate A**

Request independent specification and quality reviews. Require explicit approval
of symmetry completeness, collective contracts, workspace receipts, and absence
of old global WPW state before continuing.

**Step 6: Commit**

```bash
git add CMakeLists.txt src/common/dg_hybrid_windowed_pw_basis.f90 tests/dg/test_dg_hybrid_windowed_pw_basis_mpi.f90 tests/dg/run_dg_hybrid_windowed_pw_basis_mpi.py
git commit -m "feat: build covariant windowed PW packets"
```

### Task 4: Implement The Local Wannier Orthogonal Complement

**Files:**
- Create: `src/common/dg_hybrid_wannier_complement.f90`
- Create: `tests/dg/test_dg_hybrid_wannier_complement_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_wannier_complement_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write the dense-oracle failing test**

Construct complex orthonormal Wanniers and overlapping windowed PWs. Compare

```text
p_perp = p - W(W^H p)
```

against dense direct projection. Test nonlocal tails deliberately above and below
the cutoff, uneven row ownership, local row permutation, full symmetry packets,
same-rank duplicate IDs, nonorthonormal Wannier rejection, and MPI-invariant
fingerprints.

**Step 2: Verify failure on MPI 1/2/4/8**

**Step 3: Implement bounded local projection**

Store only packet-to-nearby-Wannier overlaps. Stream remote Wannier rows or use a
bounded owner schedule. Add a full-cell diagnostic mode that measures the omitted
tail, but never persist `N_W x N_P` coefficients.

**Step 4: Add an old-formulation oracle**

For the fixture only, reconstruct the historical dense formula and require equal
projected subspace/projector within tolerance. This is not production code.

**Step 5: Verify and commit**

Commit: `feat: project windowed PW packets out of Wannier space`.

### Task 5: Assemble And Rank-Reveal Sparse S_PP

**Files:**
- Create: `src/common/dg_hybrid_sparse_metric.f90`
- Create: `tests/dg/test_dg_hybrid_sparse_metric_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_sparse_metric_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing metric tests**

Verify sparse apply equals dense `S_PP`, Hermiticity, positive numerical rank,
packet-atomic pivoting, symmetry covariance, and condition receipt. Include exact
linear dependence, near-threshold complete clusters, indefinite corruption,
missing reverse edge, metadata disagreement, and allocation/extent adverse cases.

**Step 2: Implement sparse graph and packet rank selection**

Keep `S_PP`; do not compute or store `S_PP^-1/2`. Return a stable active packet
catalog and block-local preconditioner metadata.

**Step 3: Verify MPI 1/2/4/8 and commit**

Commit: `feat: retain sparse metric for projected PW packets`.

### Task 6: Assemble S/H/Z With A Full-Cell Oracle

**Files:**
- Create: `src/common/dg_hybrid_sparse_operators.f90`
- Create: `src/gs/dc/dg_hybrid_full_cell_operator_adapter.f90`
- Create: `tests/dg/test_dg_hybrid_sparse_operators_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_sparse_operators_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing operator tests**

Build a small complex WW/WP/PP basis. Compare sparse blocks with direct full-cell
quadrature for metric, `hpsi` Hamiltonian, and all three position components.
Require Hermiticity, expected W-P metric zero, nonzero W-P position transition,
symmetry covariance, and gauge covariance. Add a nonlocal pseudopotential fixture
or a callback seam proving the total `hpsi` path is used.

**Step 2: Implement the bounded adapters**

Reuse `project_dg_full_cell_hamiltonian_tiles` patterns and the existing `hpsi`
callback. Materialize one basis tile at a time. Store only sparse owned/neighbor
blocks after oracle comparison.

**Step 3: Add provenance and memory receipts**

Bind selection, window, packet, complement, metric, Hamiltonian, and position
fingerprints. Report persistent and transient peaks separately.

**Step 4: Verify**

Run hybrid MPI suites plus existing full-cell, W90, construction, and route suites.

**Step 5: Review Gate B**

Obtain independent reviews of operator orientation, periodic position convention,
nonlocal Hamiltonian coverage, sparse-vs-full oracle, collective safety, and peak
memory accounting. Resolve every Critical/Important item.

**Step 6: Commit**

Commit: `feat: assemble sparse hybrid metric Hamiltonian and position`.

### Task 7: Add A Matrix-Free Generalized Metric Solve

**Files:**
- Create: `src/rt/dg/rt_dg_hybrid_metric_solver.f90`
- Create: `tests/dg/test_rt_dg_hybrid_metric_solver_mpi.f90`
- Create: `tests/dg/run_rt_dg_hybrid_metric_solver_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing solver fixtures**

Compare matrix-free `S x`, preconditioned `S y=b`, and residuals against dense
LAPACK on small systems. Cover multiple RHS, complex coefficients, uneven packet
ownership, near-cutoff accepted metric, singular rejected metric, iteration cap,
rank disagreement, and allocation failure cleanup.

**Step 2: Implement the minimal solver**

Use sparse block apply and a block-local factor/preconditioner. Do not copy the
historical global `a_owned_w_global_p` solver. Make convergence tolerance relative
to the metric condition receipt and report iterations/residual.

**Step 3: Verify and commit**

Commit: `feat: solve the sparse hybrid metric matrix free`.

### Task 8: Add The Generalized Length-Gauge Propagator

**Files:**
- Create: `src/rt/dg/rt_dg_hybrid_length_gauge.f90`
- Create: `tests/dg/test_rt_dg_hybrid_length_gauge_mpi.f90`
- Create: `tests/dg/run_rt_dg_hybrid_length_gauge_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing propagation tests**

For a dense oracle system, compare one and many steps of

```text
i S dc/dt = (H0 + E.Z) c
```

with dense generalized propagation. Verify zero-field metric norm and energy,
field-driven polarization, WW/WP/PP cross terms, complex-gauge covariance,
permutation invariance, and branch-continuous polarization. Include a test that
fails if propagation and observable use different `Z`.

**Step 2: Implement the propagator and observables**

Keep the basis fixed. Apply sparse `H0`, `Z`, and `S` through one shared catalog.
Use the metric solver from Task 7. Start with a correctness-first Krylov/exponential
action; optimize only after equivalence tests pass.

**Step 3: Verify MPI 1/2/4/8 and commit**

Commit: `feat: propagate hybrid coefficients in the length gauge`.

### Task 9: Add Checkpoint And Restart Provenance

**Files:**
- Create: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Create: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`
- Create: `tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing round-trip tests**

Checkpoint catalog, sparse operators, coefficient ownership, metric rank, and the
complete fingerprint chain. Restart under a different MPI decomposition and
compare coefficients/observables. Add stale-window, changed selection, incomplete
packet, and corrupt fingerprint REDs.

**Step 2: Implement versioned checkpointing**

Do not reuse deleted WPW checkpoint formats. Keep the format explicitly versioned
and fail on incompatible provenance.

**Step 3: Review Gate C**

Review generalized propagation, norm/observable algebra, checkpoint MPI
independence, and lack of hidden dense/global state. Resolve findings.

**Step 4: Verify and commit**

Commit: `feat: checkpoint the windowed PW hybrid RT state`.

### Task 10: Connect A Diagnostic-Only GS Production Route

> **Superseded on 2026-08-21:** The approved ground-state route is the Si64
> self-consistent design in
> `2026-08-21-si64-hybrid-self-consistent-ground-state-design.md`.  A Si8
> one-shot generalized eigensolve remains useful only as an algebraic fixture;
> it must not publish an RT initial state.

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/check_dg_hybrid_windowed_pw_route.py`
- Create: `tests/dg/run_dg_hybrid_windowed_pw_gs_mpi.py`

**Step 1: Write the failing route/integration test**

Require ordering:

```text
Wannier completion -> block selection -> windows/G stars -> complement
-> sparse S/H/Z -> oracle comparison -> checkpoint publication
```

Forbid RT selection changes and forbidden global arrays.

**Step 2: Add an explicit diagnostic flag**

The first production connection builds and validates the hybrid basis but does
not replace the existing RT route. Default remains off. Add collective input
validation and clear receipts.

**Step 3: Retain Si8 only as a full-cell algebraic oracle**

Use MPI=8, OMP=1. Compare static generalized spectrum and operator receipts with
the current full-cell Wannier/reference calculation. Do not treat this as a
material acceptance test and do not publish an RT checkpoint.

**Step 4: Verify and commit**

Commit: `feat: connect diagnostic windowed PW hybrid construction`.

### Task 11: Connect Hybrid RT Behind An Explicit Flag

> **Additional prerequisite:** RT connection is blocked until the Si64 hybrid
> SCF has converged and the checkpoint stores the complete distributed occupied
> manifold plus occupations.  The former single-vector checkpoint is only a
> propagation primitive fixture.

**Files:**
- Modify: `src/rt/initialization_rt.f90`
- Modify: `src/rt/main_tddft.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/check_dg_hybrid_windowed_pw_route.py`
- Create: `tests/dg/run_dg_hybrid_windowed_pw_rt_mpi.py`

**Step 1: Write a failing end-to-end short-RT test**

Compare hybrid RT to a dense generalized reference for zero field and a weak
pulse. Verify norm, energy, polarization, current, and restart. Assert the basis
catalog fingerprint is unchanged at every step.

**Step 2: Connect the explicit route**

Consume the published hybrid checkpoint. Do not rebuild, reclassify, or
reorthogonalize the basis during RT.

**Step 3: Run MPI 1/2/4/8 and commit**

Commit: `feat: connect windowed PW hybrid real-time propagation`.

### Task 12: Material And Scaling Acceptance

**Files:**
- Create: `tests/dg/run_dg_hybrid_windowed_pw_acceptance.py`
- Create: `docs/verification/windowed-pw-hybrid-acceptance.md`
- Modify: `tests/dg/check_dg_hybrid_windowed_pw_route.py`

**Step 1: Add receipt parsing before long runs**

Parse accepted/rejected block counts, PW packet count, maximum graph degree,
`S_PP` condition estimate, persistent/transient bytes, communicated bytes, metric
iterations, and operator time. Abort the test if receipts are absent.

**Step 2: Run staged systems**

Use MPI=8 and OMP=1 unless the user explicitly changes it:

1. Si8 correctness smoke;
2. Si64 through hybrid construction and short RT;
3. one oxide repeated-shell case;
4. one water/solution case.

Monitor process RSS as an external diagnostic, but base acceptance on internal
owned-memory receipts and numerical checks.

**Step 3: Assess locality honestly**

Report packet count and neighbor degree versus system size. Claim linear scaling
only if both remain bounded. Otherwise keep the full-cell oracle and document the
observed scaling without adding fallback complexity.

**Step 4: Review Gate D — final review**

Request full specification and quality reviews covering physics, symmetry,
length-gauge consistency, MPI safety, memory, communication, restart, and tests.
Resolve all Critical/Important findings.

**Step 5: Run final verification**

```bash
python3 tests/dg/check_dg_hybrid_windowed_pw_route.py
python3 tests/dg/run_dg_hybrid_wannier_selection_mpi.py
python3 tests/dg/run_dg_hybrid_windowed_pw_basis_mpi.py
python3 tests/dg/run_dg_hybrid_wannier_complement_mpi.py
python3 tests/dg/run_dg_hybrid_sparse_metric_mpi.py
python3 tests/dg/run_dg_hybrid_sparse_operators_mpi.py
python3 tests/dg/run_rt_dg_hybrid_metric_solver_mpi.py
python3 tests/dg/run_rt_dg_hybrid_length_gauge_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_dg_hybrid_windowed_pw_gs_mpi.py
python3 tests/dg/run_dg_hybrid_windowed_pw_rt_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
git diff --check
```

Expected: every focused runner passes MPI 1/2/4/8; route and diff checks pass.

**Step 6: Commit**

```bash
git add docs/verification/windowed-pw-hybrid-acceptance.md tests/dg/run_dg_hybrid_windowed_pw_acceptance.py tests/dg/check_dg_hybrid_windowed_pw_route.py
git commit -m "test: validate windowed PW hybrid acceptance"
```

## Completion Criteria

The feature is complete only when:

- symmetry blocks and G stars are handled atomically;
- W-P metric overlap is below tolerance without a persistent global projection
  matrix;
- sparse `S/H/Z` agree with full-cell oracles;
- length-gauge propagation and observables use the identical `Z`;
- restart is MPI-decomposition independent;
- no basis adaptation occurs during RT;
- internal peak memory and communication receipts are bounded and reviewed;
- Si64 reaches and passes the post-Wannier hybrid construction and short-RT gates;
- all four independent review gates have no Critical or Important findings.
