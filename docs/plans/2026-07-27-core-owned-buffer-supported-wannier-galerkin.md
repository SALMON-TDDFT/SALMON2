# Core-Centered Buffer-Supported Wannier Galerkin Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a default-off overlapping Wannier Galerkin ground-state and coefficient-only RT route whose center ownership is fragment-core based, whose retained orbital tails extend through the buffer, and whose physical quadrature counts every core point exactly once.

**Architecture:** Conventional fragment DC supplies the buffer-periodic occupied-plus-empty candidate window. Center-owned Wannier functions retain their buffer tails, while overlap, Hamiltonian, density, and observable matrices are assembled over the unique global core partition. The fixed generally nonorthogonal basis solves \(H[n]C=SC\varepsilon\), and accepted real time solves \(iS\dot C=H(t)C\).

**Tech Stack:** Fortran 2008, MPI, SALMON DC/Poisson/XC/pseudopotential infrastructure, LAPACK for focused fixtures, distributed orthogonalized coefficient iteration for production, Python source contracts, CMake.

---

## Non-negotiable constraints

- Work only in:
  `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`
- Keep the parent worktree read-only. Do not copy its uncommitted implementation.
- Preserve commits `9438312`, `08b4fab`, `8bee86d`, and `05e63f7` as historical default-off direct-DG work.
- Inspect every uncommitted Task 5 diagnostic change. Retain only code needed by the new route or by reproducible failure evidence.
- Use a new default-off route flag. Do not silently reinterpret the existing direct-DG flag.
- The new route must not call direct real-space DG continuation, LCFO, EigenExa, WPW, normal checkpoint publication, or conventional RT.
- With the new route disabled, normal DC LCFO plus EigenExa behavior must remain unchanged.
- Every task requires:
  1. recorded RED;
  2. minimal implementation;
  3. focused verification;
  4. specification review;
  5. code-quality review;
  6. resolution of every Critical and Important finding;
  7. `git diff --check`;
  8. clean-first parent-prerequisite overlay build without `rsync --delete`;
  9. scoped commit.
- Do not start RT tasks unless the Si64 overlapping-Wannier ground-state gate passes.

## Overlay build protocol

For every task:

1. Create a fresh temporary source and build root with `mktemp -d`.
2. Copy the committed branch source into the source root.
3. Overlay only the enumerated parent prerequisite files needed by the
   current branch; never overlay the parent's entire dirty tree.
4. Overlay the task files from this worktree.
5. Configure with the verified MPI/EigenExa-capable toolchain, but keep
   `yn_eigenexa='n'` on every new-route run.
6. Run:

```bash
cmake --build <fresh-build-directory> -j4
```

Expected: `[100%] Built target salmon`.

Record the source manifest and build command in the task evidence.

### Task 1: Quarantine direct-DG diagnostics and define the new route

**Files:**

- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/scf_iteration_dft.f90`
- Modify: `src/gs/dc/dcdft.f90`
- Create: `tests/dg/check_dg_overlapping_wannier_route.py`
- Create: `docs/plans/2026-07-27-direct-dg-task5-diagnostic-record.md`

**Step 1: Inventory the uncommitted Task 5 diff**

Run:

```bash
git status --short
git diff -- src/common/dg_dc_direct_sipg.f90 \
  src/gs/conjugate_gradient.f90 \
  src/gs/dc/dcdft.f90 \
  src/gs/dc/dg_dc_buffer_core_faces.f90 \
  src/gs/dc/dg_dc_local_basis_ground_state.f90 \
  src/gs/scf_iteration.f90 \
  src/gs/subspace_diagonalization.f90 \
  src/io/inputoutput.f90 \
  src/io/salmon_global.f90
```

Record:

- the original full-face double-counting failure;
- the local-buffer reference subtraction result;
- the stable continuation through approximately
  \(\lambda=2.6\times10^{-2}\);
- the increasing Hermiticity defect;
- projector residual/escape failures;
- exact log paths and hashes;
- that no accepted Task 5 checkpoint was published.

Do not claim that Task 5 passed.

**Step 2: Write the RED route contract**

The checker must require a new flag such as:

```fortran
yn_dg_dc_overlapping_wannier = 'n'
```

and prove:

- default is `n`;
- only Gamma, real, non-SOI, LDA, DC input is accepted;
- enabling the new route rejects LCFO, EigenExa, WPW, direct-DG,
  conventional RT, and normal checkpoint settings;
- disabling the route leaves the normal DC call graph unchanged.

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because the route flag and dispatch do not exist.

**Step 3: Implement minimal route isolation**

Add the flag, input validation, and an explicit route dispatcher that stops
with:

```text
overlapping Wannier route: construction not implemented
```

before any forbidden stage can execute.

Do not add basis or solver behavior in this task.

**Step 4: Focused verification**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS.

Run one normal DC fixture with the new flag absent and confirm its LCFO and
EigenExa markers are unchanged.

**Step 5: Reviews, overlay, and commit**

Review route isolation and diagnostic accuracy. Resolve all
Critical/Important findings, run the fresh overlay build, then commit only
Task 1 files:

```bash
git commit -m "feat(dg): isolate overlapping Wannier route"
```

### Task 2: Define center-owned buffer-supported basis metadata

**Files:**

- Create: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_overlapping_wannier_metadata_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_metadata_mpi.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write RED metadata tests**

Define test cases with two and four fragments. Require:

- one authoritative center owner per Wannier;
- unique physical core ownership;
- buffer-supported value and gradient arrays;
- global physical IDs for every retained tail point;
- candidate, target, and retained ranks;
- basis generation and geometry fingerprints;
- rejection of duplicate center owners;
- rejection of duplicate or missing core ownership;
- rejection of stale or incomplete tail maps;
- acceptance of ranks with zero owned core points when they are valid
  buffer-only MPI shards.

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_metadata_mpi.py
```

Expected: FAIL because the metadata module does not exist.

**Step 2: Implement the metadata**

Create types equivalent to:

```fortran
type :: s_dg_wannier_tail
  integer :: center_fragment=0
  integer :: center_owner_rank=-1
  integer :: local_index=0
  integer :: generation=0
  integer(8) :: geometry_fingerprint=0_8
  integer(8),allocatable :: physical_grid_ids(:)
  real(8),allocatable :: value(:)
  real(8),allocatable :: gradient(:,:)
end type

type :: s_dg_overlapping_wannier_basis
  integer :: candidate_rank=0
  integer :: target_rank=0
  integer :: retained_rank=0
  integer :: generation=0
  type(s_dg_wannier_tail),allocatable :: tail(:)
end type
```

Use checked 64-bit products before allocations and collective validation
before any MPI payload exchange.

**Step 3: Focused verification**

Run the fixture on 1, 2, 4, and 8 ranks. Expected: PASS with identical
fingerprints and ownership counts.

Also run Task 1 contracts and `git diff --check`.

**Step 4: Reviews, overlay, and commit**

Review physical-ID canonicalization, zero-core ranks, overflow checks, and
collective failure behavior. Resolve all Critical/Important findings,
overlay-build, then commit:

```bash
git commit -m "feat(dg): add overlapping Wannier basis metadata"
```

### Task 3: Assemble the unique-core overlap metric

**Files:**

- Create: `src/gs/dc/dg_overlapping_wannier_metric.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_overlapping_wannier_metric_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_metric_mpi.py`

**Step 1: Write the RED algebra test**

For overlapping tails \(W\), compare the distributed result with:

\[
S=W^\dagger M_{\rm core}W.
\]

Require:

- every physical core point contributes once;
- all tails covering that point contribute;
- off-fragment overlap blocks are retained;
- \(S=S^\dagger\);
- sign, permutation, and unitary candidate-window invariance;
- positive-metric rank revelation;
- rejection of duplicate core quadrature;
- rejection of missing owner pairs;
- deterministic results on 1, 2, 4, and 8 ranks.

Run the new runner. Expected: FAIL.

**Step 2: Implement owner-pair assembly**

At each uniquely owned core point:

1. gather every tail whose physical-ID list covers the point;
2. form the local outer product;
3. accumulate into canonical global Wannier indices;
4. reduce once over the total route communicator.

Do not integrate buffer points independently.

Return the retained eigenspace of \(S\), its minimum eigenvalue, condition
number, and rejected null-space rank. Do not call EigenExa.

**Step 3: Focused verification**

Run the MPI fixture under bounds/FPE checking and compare bitwise-stable
ownership counts plus tolerance-stable matrices. Expected: PASS.

**Step 4: Reviews, overlay, and commit**

Review the global index map, Hermitian accumulation, metric threshold, and
memory scaling. Resolve all Critical/Important findings, overlay-build, then
commit:

```bash
git commit -m "feat(dg): assemble unique-core Wannier metric"
```

### Task 4: Construct and localize buffer-supported Wannier functions

**Files:**

- Create: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90`
- Create: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

**Step 1: Write RED construction tests**

Use candidate orbitals with:

- sign changes;
- permutations;
- a degenerate rotated subspace;
- occupied plus empty candidates;
- tails crossing one and two neighboring cores.

Require invariant retained space, deterministic center ownership, requested
target rank, occupied inclusion, complete value/gradient tails, and no
abrupt core truncation.

Add buffer-boundary value and gradient diagnostics.

Run the new runner. Expected: FAIL.

**Step 2: Implement construction**

Consume the conventional fragment DC candidate window before LCFO or
EigenExa. Localize within each auxiliary periodic fragment box, select
center-owned functions, and retain the complete configured buffer tail.

Apply the localization transform to values and gradients. Map every tail
sample to a physical global ID. Then assemble the global core metric from
Task 3 and remove only its numerical null space.

**Step 3: Add construction gates**

Reject:

- candidate or target rank loss;
- occupied inclusion above tolerance;
- nonpositive retained metric;
- tail value or gradient above the buffer-boundary tolerance;
- inconsistent center ownership or transform fingerprints.

**Step 4: Focused verification**

Run construction and metric fixtures on 1, 2, 4, and 8 ranks, route
contracts, bounds/FPE checks, and `git diff --check`.

**Step 5: Reviews, overlay, and commit**

Review gauge invariance, localization convention, occupied inclusion,
tail retention, and absence of LCFO/EigenExa/WPW calls. Resolve all
Critical/Important findings, overlay-build, then commit:

```bash
git commit -m "feat(dg): construct buffer-supported Wannier basis"
```

### Task 5: Assemble symmetric weak kinetic and local-potential matrices

**Files:**

- Create: `src/gs/dc/dg_overlapping_wannier_operators.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_overlapping_wannier_operators_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_operators_mpi.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write RED operator tests**

Compare distributed matrices with direct unique-core references:

\[
T_{ij}
=
\frac12\sum_K\int_{\Omega_K}\nabla w_i^*\cdot\nabla w_j,
\]

\[
V_{ij}
=
\sum_K\int_{\Omega_K}w_i^*Vw_j.
\]

Require:

- off-fragment tail blocks;
- \(T=T^\dagger\);
- \(V=V^\dagger\);
- constant and plane-wave kinetic references;
- sign, permutation, and retained-space rotation covariance;
- no independent direct-SIPG face addition;
- no buffer-volume double counting;
- equivalence across MPI decompositions.

Run the new runner. Expected: FAIL.

**Step 2: Implement one weak bilinear form**

Use unique-core quadrature for both integration and ownership. Use buffer
tails only to supply values and gradients at the owned core point.

If a residual discontinuity requires a face term, introduce it only through
one symmetric weak-form routine whose canonical face owner accumulates both
conjugate blocks. The default accepted continuous-tail path has no separate
face correction.

**Step 3: Project the existing local potential**

Reuse the conventional DC `vlocal` fields on unique cores. Do not recompute
ionic, Hartree, or XC physics in the new module.

**Step 4: Focused verification**

Run operator, metric, construction, and route tests; bounds/FPE checks; and
`git diff --check`. Expected: PASS.

**Step 5: Reviews, overlay, and commit**

Review integration by parts, boundary handling, double-count prevention,
Hermiticity, and stencil/tail coverage. Resolve every Critical/Important
finding, overlay-build, then commit:

```bash
git commit -m "feat(dg): assemble weak Wannier Hamiltonian"
```

### Task 6: Add nonlocal pseudopotential and observable matrices

**Files:**

- Create: `src/gs/dc/dg_overlapping_wannier_nonlocal.f90`
- Create: `src/gs/dc/dg_overlapping_wannier_observables.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_overlapping_wannier_nonlocal_mpi.f90`
- Create: `tests/dg/test_dg_overlapping_wannier_observables_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_physical_matrices_mpi.py`

**Step 1: Write RED nonlocal tests**

Use overlapping atom-projector supports. Require unique atom/projector
ownership, complete tail-projector overlaps, and
\(V_{\rm NL}=V_{\rm NL}^\dagger\). Reject duplicate atom contributions and
incomplete tails.

**Step 2: Write RED observable tests**

Require:

- position matrices with explicit cell/origin convention;
- canonical momentum diagnostic;
- physical velocity from the complete discrete commutator
  \(v=i[H,r]\);
- explicit nonlocal commutator contribution;
- gauge covariance;
- Hermiticity of position and velocity;
- the appropriate anti-Hermiticity relation before multiplying derivatives
  by \(-i\);
- agreement with direct small-system references.

Run the runner. Expected: FAIL.

**Step 3: Implement unique-owner physical matrices**

Accumulate every atom/projector once. Use the same core quadrature and tail
owner map as \(S\), \(T\), and \(V\). Do not use canonical momentum as an
automatic replacement for physical velocity.

**Step 4: Focused verification**

Run all Task 6 fixtures on 1, 2, 4, and 8 ranks, plus previous focused tests,
bounds/FPE checking, and `git diff --check`.

**Step 5: Reviews, overlay, and commit**

Review atom ownership, nonlocal support, origin conventions, commutator
signs, and gauge covariance. Resolve all Critical/Important findings,
overlay-build, then commit:

```bash
git commit -m "feat(dg): add Wannier physical matrices"
```

### Task 7: Solve the nonorthogonal coefficient problem and reconstruct density

**Files:**

- Create: `src/gs/dc/dg_overlapping_wannier_solver.f90`
- Create: `src/gs/dc/dg_overlapping_wannier_density.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_overlapping_wannier_solver_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_solver_mpi.py`

**Step 1: Write RED generalized-solver tests**

Use Hermitian \(H\) and positive nonidentity \(S\). Require:

\[
HC=SC\varepsilon,
\qquad
C^\dagger SC=I.
\]

Check degenerate rotations, occupied and empty targets, MPI distribution,
residual convergence, deterministic failure on indefinite/ill-conditioned
metrics, and a dense LAPACK reference used only by the fixture.

**Step 2: Write RED density tests**

Require:

\[
\rho(\mathbf r)
=
\sum_{ij}w_i(\mathbf r)P_{ij}w_j^*(\mathbf r)
\]

on unique cores, including all covering tails, and:

\[
\operatorname{Tr}(PS)
=
\int_\Omega\rho
=N_e.
\]

Reject duplicate core writes and missing tail pairs.

Run the runner. Expected: FAIL.

**Step 3: Implement the distributed solver**

Use metric orthogonalization and the existing coefficient-CG primitives
where applicable. Hold \(H\), \(S\), and the basis generation fixed for one
inner solve. Do not call EigenExa.

**Step 4: Implement density reconstruction**

Write density only on the rank's uniquely owned core points. Reduce to the
existing total-DC density ownership without publishing buffer density.

**Step 5: Focused verification**

Run solver and density fixtures on 1, 2, 4, and 8 ranks; previous focused
tests; bounds/FPE checks; and `git diff --check`.

**Step 6: Reviews, overlay, and commit**

Review metric orthogonalization, residual definition, occupation handling,
tail coverage, and electron-count identities. Resolve all
Critical/Important findings, overlay-build, then commit:

```bash
git commit -m "feat(dg): solve overlapping Wannier coefficients"
```

### Task 8: Add self-consistent Wannier ground state and route checkpoint

**Files:**

- Create: `src/gs/dc/dg_overlapping_wannier_scf.f90`
- Create: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90`
- Create: `tests/dg/test_dg_overlapping_wannier_scf_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_scf_mpi.py`
- Create: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_checkpoint_mpi.py`

**Step 1: Write RED SCF tests**

Use a two-fragment density-dependent toy model. Require:

- fixed basis and \(S\) during inner solves;
- \(H[n]\) rebuild only at outer SCF boundaries;
- core-only density;
- seeded mixing history;
- unmixed fixed-point verification;
- atomic rollback of coefficients, density, potentials, and fingerprints;
- explicit failure without fallback.

**Step 2: Write RED checkpoint tests**

Require a route-specific manifest and rank shards containing:

- basis and geometry generations;
- center ownership and tail physical IDs;
- \(S\), coefficients, occupations, and density;
- operator fingerprints and acceptance metrics.

Verify write/read/write equivalence and reject stale, partial, normal-DC, or
direct-DG checkpoints.

**Step 3: Implement the outer SCF**

Reuse existing DC Hartree, XC, local-potential, and mixing implementations.
Do not call direct real-space DG, LCFO, EigenExa, or WPW.

Publish the new checkpoint only after every focused ground-state gate passes.

Before publication, connect a concrete `main_dft`-owned production adapter.
Do not use an optional global procedure registration. The adapter must:

- consume the real DC buffer-periodic candidate orbitals and unique-core
  physical state;
- pass the periodic box IDs, symmetry image map, localization phases, values,
  and gradients through construction;
- recompute symmetry closure from those authoritative mapped values and
  gradients and bind it to the basis fingerprint;
- assemble the weak unique-core operators, execute SCF, and publish/reuse the
  route-specific checkpoint within one rollback-capable transaction.

Add RED route-level tests proving that an accepted transaction publishes,
SCF or unmixed failure preserves the prior manifest, one-sided periodic-tail
corruption fails collectively, and the new route cannot reach normal shutdown
checkpoint publication.

**Step 4: Focused verification**

Run both fixtures on 1, 2, 4, and 8 ranks, all previous focused tests,
bounds/FPE checks, normal-DC regression, and `git diff --check`.

**Step 5: Reviews, overlay, and commit**

Review SCF ordering, unmixed gate, rollback completeness, checkpoint
ownership, and route isolation. Resolve all Critical/Important findings,
overlay-build, then commit:

```bash
git commit -m "feat(dg): converge overlapping Wannier ground state"
```

### Task 9: Pass the Si64 overlapping-Wannier ground-state gate

**Files:**

- Create: `tests/dg/check_si64_overlapping_wannier_gate.py`
- Create: DG-only Si64 input fixtures under the existing sample/test-data convention
- Create: `docs/plans/2026-07-27-si64-overlapping-wannier-results.md`

**Step 1: Record the unchanged normal reference**

Use Gamma, non-SOI, PZ, `nstate_frag=400`, `mixrate=0.2`,
`threshold=1e-9`, 8 MPI, and 1 OpenMP with the new route disabled.

Confirm the known SCF trace and normal LCFO plus EigenExa markers.

**Step 2: Write the RED evidence checker**

The checker must derive acceptance metrics from raw logs and checkpoint
contents, not trust hand-authored Boolean fields.

Require the exact Cartesian matrix:

- two equivalent fragment decompositions;
- two buffer widths;
- two candidate/target Wannier windows.

For every row parse and validate runtime, memory, ranks, metric spectrum,
centers, spreads, tail value/gradient norms, \(S/H/T/V_{\rm NL}/v\)
defects, density and coefficient residuals, charge, energy, and checkpoint
round trip.

Run:

```bash
python3 tests/dg/check_si64_overlapping_wannier_gate.py <fresh-result-root>
```

Expected: FAIL because no accepted matrix exists.

**Step 3: Run the matrix**

Use fresh directories. No run may contain direct-DG, LCFO, EigenExa, WPW,
normal checkpoint, or conventional RT markers.

Record the complete raw evidence and hashes.

**Step 4: Enforce the gate**

At least one configuration must pass:

- occupied inclusion;
- positive well-conditioned \(S\);
- buffer tail value and gradient gates;
- buffer/window convergence of \(S,H,T,v,\rho\);
- all operator Hermiticity gates;
- generalized eigen residual and \(S\)-orthonormality;
- unmixed density fixed point;
- \(\operatorname{Tr}(PS)\), integrated density, and 256 electrons;
- route checkpoint round trip.

If none passes:

- record the failure evidence;
- mark Task 9 blocked;
- do not start Task 10;
- do not publish to forbidden routes.

**Step 5: Reviews, overlay, and commit**

If a configuration passes, perform specification and quality reviews
against the raw artifacts, resolve every Critical/Important finding, run a
fresh overlay build, then commit:

```bash
git commit -m "test(dg): pass Si64 overlapping Wannier gate"
```

### Task 10: Propagate coefficients in the fixed nonorthogonal basis

**Prerequisite:** Task 9 passed. Otherwise do not start.

**Files:**

- Create: `src/rt/dg/rt_dg_overlapping_wannier.f90`
- Modify: `src/rt/CMakeLists.txt`
- Modify: `src/rt/main_tddft.f90`
- Create: `tests/dg/test_rt_dg_overlapping_wannier_mpi.f90`
- Create: `tests/dg/run_rt_dg_overlapping_wannier_mpi.py`

**Step 1: Write RED metric-propagation tests**

Require propagation of:

\[
iS\dot C=H(t)C
\]

with fixed basis, fixed \(S\), and route checkpoint fingerprints.

Check:

- \(C^\dagger SC\) conservation;
- field-free energy conservation;
- phase and retained-space rotation covariance;
- time-dependent local field coupling through validated position/velocity
  matrices;
- short-run and restart-split equivalence;
- explicit rejection of stale basis, tail escape, or operator mismatch.

Run the runner. Expected: FAIL.

**Step 2: Implement coefficient-only RT**

Load only the accepted Task 9 checkpoint. Hold Wannier values, gradients,
ownership, and \(S\) fixed. Update the coefficient Hamiltonian using the
validated observable matrices.

Do not call conventional orbital RT and do not update the basis
automatically.

**Step 3: Focused verification**

Run on 1, 2, 4, and 8 ranks with bounds/FPE checking, plus every focused GS
test and `git diff --check`.

**Step 4: Reviews, overlay, and commit**

Review metric propagation, field coupling, norm/energy conservation,
restart identity, and forbidden-route isolation. Resolve all
Critical/Important findings, overlay-build, then commit:

```bash
git commit -m "feat(dg): propagate overlapping Wannier coefficients"
```

## Final verification

After Task 10, and only if Task 9 passed, run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_dg_overlapping_wannier_metadata_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_metric_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_operators_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_physical_matrices_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_solver_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_scf_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_checkpoint_mpi.py
python3 tests/dg/check_si64_overlapping_wannier_gate.py <fresh-result-root>
python3 tests/dg/run_rt_dg_overlapping_wannier_mpi.py
git diff --check
```

Then run:

- a clean-first parent-prerequisite overlay build;
- the unchanged normal DC Si64 LCFO plus EigenExa reference;
- the accepted overlapping-Wannier Si64 ground state;
- short and restart-split coefficient RT.

Perform final specification and quality reviews. Resolve every
Critical/Important finding before claiming completion.
