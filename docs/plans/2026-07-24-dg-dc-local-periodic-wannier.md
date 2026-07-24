# DG-DC Local-Periodic Wannier Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a default-off Si64/LDA/Gamma/non-SOI route that hands off a loosely converged DC density to a self-consistent local-periodic real-space SIPG ground-state solve, then constructs a manifest-complete core-owned local Wannier basis.

**Architecture:** Reuse the conventional DC fragment orbitals rather than
materializing global nodal orbitals.  Orthonormalize each fragment basis in
the uniquely owned core metric, apply that same transform to its complete
periodic `core + buffer` storage, and evaluate the existing DC
kinetic/local/nonlocal Hamiltonian on the full transformed buffer basis.
Store only distributed global-band coefficient rows.  Add the complete
Hermitian SIPG coupling on the six signed physical faces, solve all occupied
and empty configured bands with the existing orthogonalized-CG/metric
Gram--Schmidt pattern, and rebuild the density from core-owned values.
Density, Hartree, XC, and local potential updates reuse the existing DC path.
No LCFO, EigenExa, global real-space wavefunction, normal checkpoint
publication, or RT promotion occurs on this route.

**Tech Stack:** Fortran 2008, MPI, SALMON DC/Poisson/XC/pseudopotential infrastructure, symmetric interior-penalty DG, LAPACK/ScaLAPACK/EigenExa, Wannier90-compatible algebra, Python source contracts, CMake.

---

## Execution boundary

- Work only in
  `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`
  on `codex/wpw-s-orthogonal-complement`.
- Expected starting HEAD: `f94404f`; expected status: clean.
- Treat `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG` as read-only.
- Do not import its uncommitted prerequisites. Full builds use a fresh `/tmp`
  overlay of parent prerequisites plus committed branch changes, without
  `--delete`.
- The committed MINRES-QLP plan is superseded and must not be executed.
- Keep the new route default `n`. Do not alter normal DC, LCFO, WPW,
  checkpoint publication, or RT behavior.
- Use `@superpowers:test-driven-development` for every behavior change,
  `@superpowers:systematic-debugging` for failures,
  `@superpowers:requesting-code-review` after every Task, and
  `@superpowers:verification-before-completion` before every commit.
- Each Task requires RED evidence, focused verification, specification review,
  code-quality review, and resolution of all Critical/Important findings.

### Authoritative Task 3--5 refinement

This section supersedes the older nodal-materialization wording inside Tasks
3--5 below.

- Task 3 hands off the accepted DC density and builds the fragment-local basis
  coefficient representation.  The core metric determines the rank and
  transform, but the transform is applied to every DC `core + buffer` point
  before the conventional DC Hamiltonian action.  Projection and density
  ownership use core points only.  The global band axis is retained in full,
  including empty states.
- The coefficient CG holds the complete Hamiltonian fixed during each inner
  solve.  Only the outer SCF rebuilds density and the reused DC
  Hartree/XC/vlocal terms.  SIPG communication covers exactly the signed six
  physical faces; edge/corner members of the DC 27-neighbor buffer halo are
  not DG faces.
- Task 4 adds adaptive continuation of the complete SIPG operator and the
  unmixed fixed-point gate in this same distributed coefficient
  representation.  It must not switch to LCFO/EigenExa or global nodal
  orbitals.
- Task 5 is the Si64/LDA/Gamma/non-SOI DG ground-state gate for this local
  representation.  Its converged coefficients and eigenvalues live in an
  explicit DG-only in-memory state; standard DC output/checkpoint publication
  remains disabled until a later task explicitly defines and verifies that
  contract.

### Task 1: Introduce a GS/RT-neutral nodal DG core

**Files:**
- Create: `src/common/dg_nodal_state.f90`
- Create: `src/common/dg_nodal_interfaces.f90`
- Modify: `src/common/CMakeLists.txt`
- Modify: `src/rt/dg/rt_dg_nodal_types.f90`
- Modify: `src/rt/dg/rt_dg_nodal_ground_state_solver.f90`
- Modify: `src/rt/dg/rt_dg_nodal_cg.f90`
- Modify: `src/rt/dg/rt_dg_nodal_density.f90`
- Create: `tests/dg/test_dg_nodal_common_mpi.f90`
- Create: `tests/dg/run_dg_nodal_common_mpi.py`
- Modify: `tests/dg/check_dg_nodal_ground_state_route.py`

**Step 1: Write RED common-state tests**

Require a common state that owns core nodal orbitals, fragment/core/buffer
geometry, face trace storage, occupations, readiness/provenance, and explicit
core ownership. Test allocation, release, malformed geometry, duplicate core
ownership, finite validation, and transactional failure on two ranks.

**Step 2: Run RED**

```bash
python3 tests/dg/run_dg_nodal_common_mpi.py
python3 tests/dg/check_dg_nodal_ground_state_route.py
```

Expected: FAIL because the common modules do not exist.

**Step 3: Implement minimal common types and interfaces**

Define:

```fortran
type, public :: s_dg_nodal_common_state
  logical :: initialized=.false., ground_state_ready=.false.
  integer :: fragment=0, nstate=0, nspin=0, halo_width=0
  integer :: core_size(3)=0, box_size(3)=0, buffer(3)=0
  integer(8) :: geometry_fingerprint=0, operator_fingerprint=0
  complex(8), allocatable :: psi_core(:,:,:,:,:)
  real(8), allocatable :: occupations(:,:)
  type(s_dg_nodal_face), allocatable :: faces(:)
end type
```

Provide abstract complete-H, density, subspace-rotation, and collective
validator callbacks. Keep RT modules as thin compatibility wrappers; no RT
behavior may change.

**Step 4: Verify**

Run the RED commands plus existing nodal GS, density, current, and Taylor
fixtures. Expected: PASS.

**Step 5: Review and commit**

Review dependency direction, ownership, cleanup, MPI agreement, and unchanged
RT semantics. Resolve all Critical/Important findings.

```bash
git commit -m "refactor(dg): introduce common nodal state core"
```

### Task 2: Implement the physical SIPG face operator

**Files:**
- Create: `src/common/dg_nodal_sipg.f90`
- Modify: `src/common/CMakeLists.txt`
- Create: `tests/dg/test_dg_nodal_sipg_mpi.f90`
- Create: `tests/dg/run_dg_nodal_sipg_mpi.py`
- Modify: `tests/dg/check_dg_nodal_ground_state_route.py`

**Step 1: Write an analytic RED fixture**

Use one-dimensional two-element polynomials embedded in the 3-D storage.
Build the dense SIPG matrix independently and compare the matrix-free action.
Test:

- consistency, symmetric, and penalty terms separately;
- Hermiticity;
- internal-face left/right cancellation;
- canonical single face ownership;
- physical-supercell periodic faces;
- auxiliary local-box periodic wrap excluded from physical coupling;
- fragment relabeling invariance;
- invalid normals, duplicate faces, stale traces, and nonfinite input.

**Step 2: Run RED**

```bash
python3 tests/dg/run_dg_nodal_sipg_mpi.py
```

Expected: FAIL because the SIPG action is absent.

**Step 3: Implement the complete face action**

For each canonical physical face implement the symmetric bilinear form:

\[
-\{\nabla_n u\}[v]-\{\nabla_n v\}[u]
+\eta h_f^{-1}[u][v].
\]

Use the existing canonical face schedule and trace exchange where possible.
Apply both element contributions in identical collective order. Treat the
consistency, symmetry, and penalty pieces as one Hermitian operator.

**Step 4: Add operator diagnostics**

Expose face Hermiticity defect, internal cancellation defect, jump norm,
penalty energy, trace epoch, and ownership counts. Fail collectively before
returning partial output.

**Step 5: Verify, review, and commit**

Run the new fixture with normal flags and `-fcheck=all -fbacktrace`, plus
canonical-face and weak-form regressions and `git diff --check`. Resolve all
Critical/Important findings.

```bash
git commit -m "feat(dg): add nodal SIPG face operator"
```

### Task 3: Add early DC handoff and full candidate materialization

**Files:**
- Create: `src/gs/dc/dg_dc_handoff.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/scf_iteration_dft.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Create: `tests/dg/test_dg_dc_handoff_mpi.f90`
- Create: `tests/dg/run_dg_dc_handoff_mpi.py`
- Create: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write RED handoff tests**

Require default-off input plumbing and a handoff state machine that accepts
only after minimum iterations, finite/correct charge, loose DC residual, and
successful fragment solves. Reject any LCFO/Wannier/window truncation before
handoff. Test thresholds `1E-2`, `1E-3`, and `1E-4`.

Require complete materialization of the configured candidate pool, core
ownership without duplicates, local-periodic `core + buffer` metadata, and
metric-rank diagnostics. Prove that post-handoff state updates can leave the
original DC coefficient span.

**Step 2: Run RED**

```bash
python3 tests/dg/run_dg_dc_handoff_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
```

Expected: FAIL for missing controls and handoff.

**Step 3: Add default-off controls**

Add:

```text
yn_dg_dc_local_periodic = 'n'
dg_dc_handoff_min_iter
dg_dc_handoff_tolerance
dg_dc_candidate_orbitals_per_atom
dg_dc_metric_rank_tolerance
```

Validate them even when disabled. For the first gate use 40 candidate
orbitals/atom and compare handoff tolerance rather than hard-coding one.

**Step 4: Implement handoff**

Stop normal DC at the first valid gate, preserve density/potentials, discard
mixing history, materialize all candidate orbitals on common nodal degrees of
freedom, validate metric rank collectively, and bind geometry/operator
fingerprints. Do not invoke LCFO or Wannier.

**Step 5: Verify, review, and commit**

Run new tests, existing DC density assembly tests, normal SCF contracts, and
`git diff --check`. Review normal-route isolation, candidate completeness,
rank policy, and cleanup.

```bash
git commit -m "feat(dg): add early DC-to-nodal handoff"
```

### Task 4: Implement adaptive self-consistent DC-to-DG continuation

**Files:**
- Create: `src/gs/dc/dg_dc_ground_state.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/common/dg_nodal_interfaces.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Create: `tests/dg/test_dg_dc_ground_state_mpi.f90`
- Create: `tests/dg/run_dg_dc_ground_state_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write RED continuation tests**

Use a two-fragment nonlinear toy potential with known fixed point. Verify:

- `lambda=0` handoff;
- intermediate SCF convergence before lambda advance;
- complete Hermitian interface scaling;
- adaptive step growth;
- rollback and step halving;
- minimum-step failure;
- occupied-projector tracking rather than state labels;
- density mixing-history reset;
- final `lambda=1`;
- unmixed fixed-point rejection and acceptance;
- stale fingerprint and nonfinite collective failure.

**Step 2: Run RED**

```bash
python3 tests/dg/run_dg_dc_ground_state_mpi.py
```

Expected: FAIL because the driver does not exist.

**Step 3: Add bounded controls**

Add controls for DC handoff, intermediate/final orbital and density
tolerances, initial/min/max lambda step, allowed residual growth, maximum
SCF/eigensolver/rollback counts, simple density mix rate, Hermiticity,
orthogonality, face-balance, and electron-count tolerances.

**Step 4: Implement the staged driver**

Implement:

```text
DC_PRECONVERGENCE
DG_CONTINUATION(lambda<1)
FULL_DG_SCF(lambda=1)
UNMIXED_FIXED_POINT
ACCEPTED / FAILED
```

At every fixed lambda, solve the complete nodal eigenproblem, reconstruct
core-owned density, update Hartree and LDA XC, mix density, and verify orbital,
density, subspace, face, and charge diagnostics. Rollback state and density
transactionally.

**Step 5: Verify, review, and commit**

Run new tests, SIPG tests, nodal complete-H/nonlocal/density regressions, and
`git diff --check`. Review continuation mathematics, rollback, MPI ordering,
mixing, and final fixed-point semantics.

```bash
git commit -m "feat(dg): add self-consistent DC-to-DG continuation"
```

### Task 5: Add DG-GS checkpoint and Si64 ground-state gate

**Files:**
- Create: `src/common/dg_ground_state_checkpoint.f90`
- Modify: `src/common/CMakeLists.txt`
- Modify: `src/gs/dc/dg_dc_ground_state.f90`
- Create: `tests/dg/test_dg_ground_state_checkpoint_mpi.f90`
- Create: `tests/dg/run_dg_ground_state_checkpoint_mpi.py`
- Modify: `docs/plans/2026-07-24-dg-dc-local-periodic-wannier.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

**Step 1: Write RED checkpoint tests**

Require transactional temporary-write/rename, schema/version, geometry,
physical-face ownership, continuation history, `lambda=1`, orbitals,
occupations, density/potentials, fingerprints, and every acceptance
diagnostic. Reject partial, stale, unaccepted, or LCFO/WPW-mislabeled data.

**Step 2: Implement and verify checkpoint**

Implement a new schema; do not modify existing LCFO/WPW checkpoint meaning.
Run round-trip, corruption, rank-disagreement, and incomplete-publication
tests.

**Step 3: Parent-prerequisite overlay build**

Record branch/HEAD/clean status and parent porcelain hash. Create a fresh
read-only-parent `/tmp` overlay containing only committed branch changes plus
temporary parent prerequisites. Configure and build:

```bash
cmake -S <overlay-source> -B <overlay-build> \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/mpifort \
  -DCMAKE_Fortran_FLAGS=-fallow-invalid-boz \
  -DUSE_MPI=ON -DUSE_SCALAPACK=ON -DUSE_EIGENEXA=ON -DUSE_WANNIER90=OFF
cmake --build <overlay-build> --target salmon -j1
```

Expected: `[100%] Built target salmon`.

**Step 4: Run Si64 ground-state matrix**

Use fresh run directories and eight MPI ranks, `OMP_NUM_THREADS=1`. Compare:

- handoff `1E-2`, `1E-3`, `1E-4`;
- at least two buffer widths;
- at least two equivalent fragment decompositions.

Record DC handoff state, metric rank, lambda history/rollbacks, orbital and
density residuals, charge, energy, orthogonality, Hermiticity, face balance,
fixed-point result, runtime, and checkpoint publication.

**Step 5: Review, document, and commit**

Review raw-log transcription and ratios. If no configuration passes all
ground-state gates, stop the plan here and do not implement Wannier Tasks 6–7.
Resolve Critical/Important findings and commit only code/tests/docs belonging
to this Task.

```bash
git commit -m "feat(dg): publish verified DG-DC ground state"
```

### Task 6: Implement atomic-orbital coverage manifests and projectability

**Prerequisite:** Task 5 has at least one accepted Si64 DG ground state.

**Files:**
- Create: `src/common/dg_orbital_coverage.f90`
- Create: `src/gs/dc/dg_atomic_projection.f90`
- Modify: `src/common/CMakeLists.txt`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_orbital_coverage.f90`
- Create: `tests/dg/run_dg_orbital_coverage.py`
- Create: `tests/dg/test_dg_atomic_projection_mpi.f90`
- Create: `tests/dg/run_dg_atomic_projection_mpi.py`

**Step 1: Write RED manifest tests**

Test s/p, s/p/d, and second-zeta shell degeneracies, species mapping,
versioning, pseudopotential fingerprinting, invalid shell/radial counts, and
deterministic tetrahedral orientation from neighbor geometry.

**Step 2: Write RED projectability tests**

Build synthetic outer eigenspaces with known missing/complete shells. Verify
projection matrix, singular spectrum, full target rank, occupied inclusion,
buffer/core ownership, and automatic outer-window expansion requests.

**Step 3: Implement minimal manifest and projection algebra**

Treat configured DC orbitals/atom as the candidate pool. Determine target
dimension from manifest shells. Remove only collective metric-null
directions. Never infer coverage solely from energy ordering or nonlocal
projector count.

**Step 4: Verify, review, and commit**

Run both new suites, generalized algebra/row-layout regressions, bounds-check
builds, and review shell semantics, rank thresholds, and scale invariance.

```bash
git commit -m "feat(dg): add atomic orbital coverage manifests"
```

### Task 7: Construct and align core-owned local Wannier functions

**Files:**
- Create: `src/gs/dc/dg_local_periodic_wannier.f90`
- Create: `src/gs/dc/dg_local_wannier_gauge.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/common/dg_ground_state_checkpoint.f90`
- Create: `tests/dg/test_dg_local_periodic_wannier_mpi.f90`
- Create: `tests/dg/run_dg_local_periodic_wannier_mpi.py`
- Create: `tests/dg/test_dg_local_wannier_gauge_mpi.f90`
- Create: `tests/dg/run_dg_local_wannier_gauge_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write RED local Wannier tests**

Verify frozen occupied-projector inclusion, manifest target rank,
disentanglement from an oversized outer window, local-periodic localization,
unique home-atom/core ownership, exclusion of buffer-only centers, retention
of Wannier tails throughout the periodic buffer, band and density
reproduction, and failure on insufficient projectability.

**Step 2: Write RED gauge-alignment tests**

For known neighboring unitary rotations, verify shared-buffer overlap,
polar/Procrustes recovery, deterministic spanning-tree alignment, fragment
relabeling invariance, and loop-closure defects. Reject rank-deficient shared
overlaps and ambiguous ownership.

**Step 3: Implement local construction**

Build manifest projections for core and buffer atoms, expand the outer window
until rank criteria pass, freeze the occupied projector, disentangle the full
manifest target, localize in each auxiliary periodic box, and retain only
core-owned centers.  Do not truncate a Wannier function at its core boundary.
Keep its translated tails throughout the auxiliary periodic buffer, recording
the center and lattice image.  The auxiliary buffer boundary is a
representation boundary, not an additional physical DG face; density is
owned once on the core and SIPG surface terms are owned once on actual
core--core interfaces.

**Step 4: Implement static gauge alignment**

Align neighboring target subspaces from shared-buffer overlaps. This is
static overlap gauge stitching only; do not add Wilson-link, field, position,
polarization, or RT code.

**Step 5: Verify, review, and commit**

Run new suites, existing Wannier generalized-algebra tests, bounds checks, and
the full overlay build. Review occupied inclusion, rank policy, locality,
ownership, loop closure, and absence of RT scope.

```bash
git commit -m "feat(dg): add core-owned local periodic Wannier basis"
```

### Task 8: Solve self-consistently in the periodic-buffer Wannier basis and run the Si64 gates

**Files:**
- Create: `src/gs/dc/dg_local_periodic_wannier_scf.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_local_periodic_wannier_scf_mpi.f90`
- Create: `tests/dg/run_dg_local_periodic_wannier_scf_mpi.py`
- Modify: `docs/plans/2026-07-24-dg-dc-local-periodic-wannier.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

**Step 1: Write RED Wannier-SCF tests**

With fixed periodic-buffer Wannier functions, verify

```text
H_W[n] = H_DC_W[n] + H_DG_W
```

is held fixed during each coefficient CG solve, including occupied and empty
bands. Reconstruct density only on uniquely owned core points, then reuse the
DC Hartree, XC, `vlocal`, and mixing updates. Require self-consistency of the
Wannier coefficients and density. Prove that auxiliary buffer wraps add no
DG interface and that increasing the periodic buffer converges tail, density,
and SIPG matrix errors.

For a complete untruncated unitary basis, verify equivalence to the accepted
Task 5 real-space/local-basis DG fixed point. For a truncated Wannier basis,
verify that a transform-only result is rejected until the Wannier-space SCF
converges.

**Step 2: Implement fixed-basis Wannier SCF**

Transform the accepted Task 5 DC-volume-plus-SIPG Hamiltonian into the
periodic-buffer Wannier basis. Keep Wannier functions fixed for the inner
SCF; solve coefficients with orthogonalized CG and metric Gram--Schmidt,
reconstruct core density, and update the existing DC potentials. Rebuild
Wannier functions only in an explicitly controlled outer loop if buffer or
subspace convergence fails. Do not invoke LCFO/EigenExa full-system
wavefunctions.

**Step 3: Run Tier 0 and Tier 1**

From an accepted Task 5 checkpoint, run:

- Tier 0: Si `3s+3p`, target 4/atom, 256 Wannier;
- Tier 1: Si `3s+3p+d`, target 9/atom, 576 Wannier;
- at least two outer-window sizes;
- at least two buffer widths.

**Step 4: Record evidence**

Record candidate/metric/outer/target ranks, projection singular values,
occupied-inclusion error, centers, spreads, ownership, neighbor overlap,
loop closure, band/density errors, pre/post-Wannier-SCF differences, Wannier
coefficient and density residuals, tail norm at the core and buffer
boundaries, buffer sensitivity of the SIPG matrix, and runtime/memory.

Select the smallest manifest tier for which requested occupied and conduction
observables, density, projectability, spreads, and retained Hamiltonian
elements are converged.

**Step 5: Final review and commit**

Review all values against raw logs. Rerun all focused tests, overlay build,
and `git diff --check`. Resolve every Critical/Important finding. Commit only
the two result documents:

```bash
git commit -m "docs(dg): record local-periodic DG-DC Wannier gates"
```

### Task 9: Decision checkpoint

Present the ground-state and Wannier evidence for user review. Do not promote
the route into normal DC, LCFO replacement, checkpoint default, publication,
or RT automatically.

If the DG ground state fails, use handoff/continuation/face diagnostics before
changing Wannier policy. If the ground state passes but Wannier fails, use
manifest projectability, outer-window rank, buffer convergence, and
gauge-alignment evidence to decide whether to enlarge the candidate pool,
manifest tier, outer window, or buffer.
