# DG Buffer-Window Projector, Wannier GS, and RT Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace direct same-state fragment face exchange with a full-window buffer-to-neighbor-core projector, converge the real-space DG ground state without LCFO or EigenExa, then build and self-consistently solve a buffer-periodic Wannier representation whose coefficients alone are propagated in real time.

**Architecture:** Conventional DC remains the owner of the initial fragment SCF, real-space volume Hamiltonian, density, Hartree, XC, and local potential. At each outer DG SCF boundary, a rank-revealing projector is built from all configured occupied and empty neighboring-core states on the physical buffer overlap; it is frozen during the inner real-space orbital CG and supplies the six-face DG traces. After the Si64 real-space gate passes, buffer-periodic Wannier orbitals are constructed, re-solved self-consistently with Wannier DG blocks, and then held fixed while only Wannier coefficients are propagated.

**Tech Stack:** Fortran 2008, MPI, OpenMP, SALMON DC/Poisson/XC/pseudopotential infrastructure, LAPACK rank-revealing factorization, symmetric interior-penalty DG, Python route contracts, CMake.

---

## Execution boundary

- Work only in
  `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`
  on `codex/wpw-s-orthogonal-complement`.
- Treat `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG` as read-only.
- Do not copy uncommitted parent-worktree implementation into this branch.
- The authoritative design is
  `docs/plans/2026-07-25-dg-buffer-window-projector-wannier-rt-design.md`.
- This plan supersedes
  `docs/plans/2026-07-25-dg-direct-dc-cg-implementation.md`.
- Preserve the useful current work for DC handoff, mixing reseeding,
  continuation snapshots, DG-only checkpointing, and route isolation.
- Remove or replace the invalid production assumption that neighboring
  fragment state `io` is the remote DG trace for local state `io`.
- Keep the route default `n`.
- With DG disabled, normal DC must retain LCFO and EigenExa behavior.
- With DG enabled, do not call LCFO, EigenExa, WPW, normal checkpoint
  publication, or conventional RT.
- Do not begin Tasks 6--8 unless Task 5 passes the Si64 real-space DG
  ground-state gate.

## Mandatory workflow for every Task

For every Task below:

1. Use `@superpowers:test-driven-development`.
2. Record the RED command and expected failure before production edits.
3. Implement the minimum behavior needed for GREEN.
4. Run the focused verification listed for that Task.
5. Run `git diff --check`.
6. Use `@superpowers:requesting-code-review` for a specification review.
7. Use `@superpowers:requesting-code-review` for a separate code-quality
   review.
8. Resolve every Critical and Important finding.
9. Re-run focused verification after review fixes.
10. Build a fresh parent-prerequisite overlay without `--delete`.
11. Use `@superpowers:verification-before-completion`.
12. Commit only the Task's reviewed files.

Use `@superpowers:systematic-debugging` before changing production code in
response to any failed or unexpected verification.

### Task 1: Implement the full-window overlap-metric projector

**Files:**

- Create: `src/common/dg_buffer_window_projector.f90`
- Modify: `src/common/CMakeLists.txt`
- Create: `tests/dg/test_dg_buffer_window_projector.f90`
- Create: `tests/dg/run_dg_buffer_window_projector.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write the failing algebra fixture**

Define a fixture with:

- a nonorthogonal restricted neighboring window \(\Phi\);
- an exactly representable buffer vector;
- a vector with a known component outside the retained window;
- independent sign changes;
- a state permutation;
- an orthogonal rotation in a degenerate two-state subspace;
- a rank-deficient overlap metric.

Require the API:

```fortran
type :: s_dg_buffer_projector_diagnostics
  integer :: configured_states=0
  integer :: retained_rank=0
  real(8) :: minimum_retained_singular_value=0d0
  real(8) :: projection_residual=huge(1d0)
  real(8) :: escape_norm=huge(1d0)
end type

subroutine build_dg_buffer_window_projector(phi_core,buffer_values,weights, &
    rank_tolerance,coefficients,reconstructed,diagnostics,ok,message)
```

Assert:

```fortran
call build_dg_buffer_window_projector(phi,buffer,weights,tol,c,reconstructed,d,ok,message)
if(.not.ok) error stop trim(message)
if(maxval(abs(reconstructed-exact)))>1d-12) error stop 'projection mismatch'
if(abs(d%projection_residual-expected_residual)>1d-12) error stop 'bad residual'
```

Run:

```bash
python3 tests/dg/run_dg_buffer_window_projector.py
```

Expected: FAIL because the module and runner do not exist.

**Step 2: Implement the projector**

Compute:

```fortran
S = matmul(transpose(phi_core*weights),phi_core)
B = matmul(transpose(phi_core*weights),buffer_values)
```

Use a symmetric rank-revealing eigendecomposition or SVD. Retain singular
values satisfying:

```fortran
sigma >= rank_tolerance*maxval(sigma)
```

Form \(S^+B\), reconstruct \(\Phi S^+B\), and report the weighted relative
residual and escape norm. Reject nonfinite inputs, negative weights, empty
windows, invalid tolerances, and collective callers that supply inconsistent
dimensions.

Do not multiply the projector by occupations.

**Step 3: Add invariance tests**

Prove that reconstructed values and diagnostics are invariant under:

- independent signs;
- a state permutation;
- an orthogonal rotation within the retained window.

Prove that an occupied-only window fails a fixture requiring a known empty
state, while the configured full window passes.

**Step 4: Focused verification**

Run:

```bash
python3 tests/dg/run_dg_buffer_window_projector.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS.

**Step 5: Review, overlay build, and commit**

Review pseudoinverse scaling, weighted algebra, allocation failures,
determinism, and the absence of occupation weighting. Resolve all
Critical/Important findings, build the fresh overlay, and commit:

```bash
git commit -m "feat(dg): add full-window buffer projector"
```

### Task 2: Map and exchange physical buffer-to-neighbor-core overlaps

**Files:**

- Create: `src/gs/dc/dg_dc_buffer_core_faces.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/dg_dc_ground_state_adapter.f90`
- Create: `tests/dg/test_dg_dc_buffer_core_faces_mpi.f90`
- Create: `tests/dg/run_dg_dc_buffer_core_faces_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write the failing MPI mapping fixture**

Use two and eight MPI ranks. Assign unique values from physical global grid
IDs and require exact equality between:

```text
fragment I positive/negative buffer points
fragment J corresponding core points
```

Test:

- all six signed faces;
- positive and negative buffer storage layout;
- physical-supercell periodic wrap;
- two buffer widths;
- tangential decomposition;
- fragment relabeling;
- invalid and duplicate physical grid IDs;
- rank disagreement before communication.

Run:

```bash
python3 tests/dg/run_dg_dc_buffer_core_faces_mpi.py
```

Expected: FAIL because the face mapper is absent.

**Step 2: Implement canonical face metadata**

Define a face object containing:

```fortran
integer :: axis,side,local_fragment,neighbor_fragment
integer :: overlap_depth,point_count,generation
integer(8),allocatable :: physical_grid_ids(:)
integer,allocatable :: local_buffer_indices(:,:),neighbor_core_indices(:,:)
```

Build exactly six signed physical faces. Reuse DC 27-neighbor data only to
resolve physical fragment neighbors; do not treat edges or corners as DG
faces.

**Step 3: Exchange the complete state window**

Exchange all locally owned orbital columns in deterministic face order.
State-distributed ranks communicate only their owned state range. Validate
global configured-state coverage collectively.

The received array must be indexed by physical overlap point and configured
state, never by an assumed local/neighbor same-state pairing.

**Step 4: Integrate the Task 1 projector**

For every signed face, build:

- neighboring restricted overlap metric;
- buffer-to-core coefficients;
- reconstructed value traces;
- reconstructed normal-derivative traces;
- rank, residual, and escape diagnostics.

Bind them to a projector generation and operator fingerprint.

**Step 5: Focused verification**

Run:

```bash
python3 tests/dg/run_dg_buffer_window_projector.py
python3 tests/dg/run_dg_dc_buffer_core_faces_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS.

**Step 6: Review, overlay build, and commit**

Review physical-point identity, six-face ownership, periodic images, MPI
ordering, all-state coverage, memory scaling, and collective failure.
Resolve all Critical/Important findings, build the fresh overlay, and commit:

```bash
git commit -m "feat(dg): map buffer-to-core face windows"
```

### Task 3: Replace same-state face exchange in the existing real-space CG

**Files:**

- Modify: `src/gs/conjugate_gradient.f90`
- Modify: `src/gs/scf_iteration.f90`
- Modify: `src/gs/scf_iteration_dft.f90`
- Modify: `src/common/dg_dc_direct_sipg.f90`
- Modify: `tests/dg/test_dg_dc_direct_cg_mpi.f90`
- Modify: `tests/dg/run_dg_dc_direct_cg_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write RED contracts for the invalid path**

Require that the production direct-DG subroutine:

- accepts a frozen buffer-window projector generation;
- uses reconstructed remote traces;
- applies the same frozen linear operator to `xk`, every `pk`, and refreshed
  `xk`;
- contains no direct `neighbor(io)` same-state trace;
- does not rebuild the projector inside `gscg_rwf`;
- keeps all `nstate_frag` occupied and empty states.

Add an MPI fixture where same-index states are sign-flipped, permuted, and
rotated. The old implementation must produce a different action; the
projected implementation must reproduce the reference action.

Run:

```bash
python3 tests/dg/run_dg_dc_direct_cg_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
```

Expected: FAIL because the production action still uses same-state neighbor
values.

**Step 2: Remove diagnostic-only same-index production code**

Delete temporary face-correlation logging and any direct same-state
neighbor trace used by the Hamiltonian. Preserve optional diagnostics only
if they are clearly labeled non-Hamiltonian and have bounded output.

**Step 3: Apply the frozen projected SIPG action**

At the outer SCF boundary, build a projector generation. Pass it read-only
through:

```fortran
solve_orbitals -> gscg_rwf -> apply_dc_direct_dg_hpsi_rwf
```

For each Hamiltonian evaluation:

```fortran
call hpsi(...)
call apply_projected_dc_dg_faces(frozen_projector,psi,dg_hpsi)
hpsi%rwf = hpsi%rwf + lambda*dg_hpsi%rwf
```

The action must remain linear for the lifetime of one orbital-CG sweep.

**Step 4: Keep existing DC ownership**

Preserve this outer ordering:

```fortran
call calc_density(...)
call calc_rho_total_dcdft(...)
call update_density_and_potential(...)
call calc_vlocal_fragment_dcdft(...)
```

Use the existing DC Gram--Schmidt and occupation path. Do not invoke a
global coefficient eigensolver.

**Step 5: Focused verification**

Run:

```bash
python3 tests/dg/run_dg_buffer_window_projector.py
python3 tests/dg/run_dg_dc_buffer_core_faces_mpi.py
python3 tests/dg/run_dg_dc_direct_cg_mpi.py
python3 tests/dg/run_dg_dc_handoff_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS.

**Step 6: Review, overlay build, and commit**

Review linearity, projector lifetime, timer balance, memory lifetime,
all-state retention, and unchanged DC density ownership. Resolve all
Critical/Important findings, build the fresh overlay, and commit:

```bash
git commit -m "fix(dg): use projected buffer traces in DC CG"
```

### Task 4: Correct continuation, rollback, and direct-state checkpointing

**Files:**

- Modify: `src/gs/dc/dg_dc_handoff.f90`
- Modify: `src/gs/scf_iteration_dft.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/common/dg_ground_state_checkpoint.f90`
- Modify: `tests/dg/test_dg_dc_handoff_mpi.f90`
- Modify: `tests/dg/test_dg_dc_direct_checkpoint_mpi.f90`
- Modify: `tests/dg/run_dg_dc_direct_checkpoint_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write RED continuation tests**

Require:

- mixing history is seeded from accepted density/Hartree/XC, never zero;
- the first iteration after a lambda increment is compared with the accepted
  pre-stage residual;
- the abnormal first residual can trigger rollback;
- projector generation is restored or invalidated on rollback;
- occupations are included in snapshots;
- all density and potential fields are restored atomically;
- lambda step and rollback limits are enforced.

Run:

```bash
python3 tests/dg/run_dg_dc_handoff_mpi.py
python3 tests/dg/run_dg_dc_direct_checkpoint_mpi.py
```

Expected: FAIL on the accepted-baseline and projector-generation contracts.

**Step 2: Implement nonzero mixing reseeding**

Use:

```fortran
call reset_dc_mixing_history_from_state(mixing,dc%rho_tot_s,dc%Vh_tot,dc%Vxc_tot)
```

Seed every input and output history slot from the accepted fields. Reset
only algorithmic counters.

**Step 3: Correct the continuation baseline**

Store:

```fortran
accepted_stage_residual
trial_stage_first_residual
```

Compare the first and subsequent trial residuals with the accepted
pre-stage baseline. Do not assign the first trial residual as the baseline.

**Step 4: Extend rollback and checkpoint state**

Persist and restore:

- fragment `core + buffer` orbitals;
- occupations;
- total and spin density;
- Hartree, XC, and local potential;
- lambda and step history;
- projector rank/residual/escape diagnostics;
- projector generation and fingerprints;
- six-face diagnostics.

Publish a DG checkpoint only after the full-scale gate. Keep normal
checkpoint publication suppressed on the DG route.

**Step 5: Focused verification**

Run:

```bash
python3 tests/dg/run_dg_dc_handoff_mpi.py
python3 tests/dg/run_dg_dc_direct_checkpoint_mpi.py
python3 tests/dg/run_dg_ground_state_checkpoint_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS.

**Step 6: Review, overlay build, and commit**

Review transactional rollback, nonzero mixing seed, generation staleness,
corruption detection, and publication isolation. Resolve all
Critical/Important findings, build the fresh overlay, and commit:

```bash
git commit -m "fix(dg): make projected continuation transactional"
```

### Task 5: Pass the Si64 real-space DG ground-state gate

**Files:**

- Modify: `docs/plans/2026-07-24-dg-dc-local-periodic-wannier.md`
- Create: `tests/dg/check_si64_dg_gate.py`
- Create or update DG-only Si64 input fixtures under the existing test-data
  convention.

**Step 1: Re-establish the normal DC reference**

Use the verified Si64 input and `atom.dat` with:

```text
Gamma
non-SOI
LDA
nstate_frag=400
mixrate=0.2
threshold=1e-9
8 MPI
1 OpenMP
```

DG must be disabled. Confirm the known normal DC trace and confirm that
LCFO+EigenExa execute.

**Step 2: Run the first projected-DG diagnostic case**

Use a fresh directory. Record before and after handoff:

- buffer-to-neighbor-core projection residual;
- rank and minimum retained singular value;
- escape norm;
- `||H_DC psi||`;
- `||lambda H_DG psi||`;
- density jump;
- energy;
- charge;
- lambda and rollback state.

If the first post-handoff iteration is abnormal, use systematic debugging
before changing any control.

**Step 3: Run the required matrix**

Run handoff tolerances:

```text
1e-2
1e-3
1e-4
1e-9 reference
```

Use at least:

- two buffer widths;
- two equivalent fragment decompositions.

Record runtime, peak memory, state count, projector time, SCF/CG counts,
continuation history, rollback history, density and orbital residuals,
energy, charge, orthogonality, Hermiticity, six-face balance, and checkpoint
state.

**Step 4: Enforce the full-scale gate**

Require at lambda one:

- unmixed density fixed point;
- orbital residual tolerance;
- orthogonality tolerance;
- projection residual tolerance;
- escape norm tolerance;
- Hermiticity tolerance;
- forward/reverse face balance;
- electron-count tolerance;
- complete six-face coverage;
- checkpoint round-trip success.

Run:

```bash
python3 tests/dg/check_si64_dg_gate.py <fresh-result-directory>
```

Expected: PASS.

**Step 5: Stop or continue**

If no Si64 configuration passes:

- record the failure evidence;
- mark Task 5 blocked;
- do not execute Tasks 6--8;
- do not publish to normal DC, LCFO, WPW, checkpoint, or RT routes.

If a configuration passes, perform both reviews, resolve all
Critical/Important findings, run the fresh overlay build, and commit:

```bash
git commit -m "test(dg): pass Si64 projected real-space gate"
```

### Task 6: Construct buffer-periodic Wannier orbitals from the accepted DG state

**Prerequisite:** Task 5 passed. Otherwise do not start.

**Files:**

- Modify: `src/gs/dc/lcfo_wannier_sawf.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf_orchestrator.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf_collective.f90`
- Modify: `src/gs/main_dft.f90`
- Create: `src/gs/dc/dg_buffer_wannier_state.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_buffer_wannier_mpi.f90`
- Create: `tests/dg/run_dg_buffer_wannier_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write RED Wannier-state tests**

Require a DG-only Wannier state containing:

- core values;
- complete retained buffer tails;
- centers and spreads;
- fragment owner and physical-periodic images;
- source `nstate_frag` window;
- gauge, geometry, and projector fingerprints.

Prove that `nstate_frag` is construction-window size and that the persisted
Wannier count is the later propagation dimension.

**Step 2: Add route-isolation RED**

Require that DG Wannier construction consumes the accepted real-space DG
checkpoint and does not call LCFO or EigenExa.

**Step 3: Implement buffer-periodic construction**

Construct Wannier orbitals on the fragment `core + buffer` periodic space.
Retain tails through the configured buffer. Validate centers, spreads, tail
norms, state-window completeness, physical image metadata, and collective
ownership.

**Step 4: Focused verification**

Run:

```bash
python3 tests/dg/run_dg_buffer_wannier_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS with no LCFO/EigenExa calls on the DG route.

**Step 5: Review, overlay build, and commit**

Resolve all Critical/Important findings and commit:

```bash
git commit -m "feat(dg): construct buffer-periodic Wannier state"
```

### Task 7: Re-solve the DG ground state self-consistently in the Wannier representation

**Prerequisite:** Tasks 5 and 6 passed.

**Files:**

- Create: `src/gs/dc/dg_wannier_ground_state.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/common/dg_ground_state_checkpoint.f90`
- Create: `tests/dg/test_dg_wannier_ground_state_mpi.f90`
- Create: `tests/dg/run_dg_wannier_ground_state_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write the failing nonlinear fixture**

Use a two-fragment Wannier toy model with density-dependent local potential.
Require:

- core-owned density reconstruction;
- reused DC Hartree/XC/vlocal update;
- updated volume and six-face DG blocks;
- iterative occupied-plus-required-empty coefficient solve;
- \(C^T S_W C=I\);
- density and coefficient fixed points;
- rollback on a failed stage;
- no LCFO or EigenExa.

**Step 2: Implement the outer Wannier SCF**

For every iteration:

```fortran
call reconstruct_wannier_core_density(...)
call update_dg_density_and_potential_with_dc(...)
call update_wannier_volume_hamiltonian(...)
call update_wannier_dg_face_blocks(...)
call solve_wannier_coefficients_iterative(...)
```

Keep the Hamiltonian fixed during each inner coefficient solve.

**Step 3: Implement the Wannier GS gate**

Require coefficient residual, density fixed point, overlap
orthonormality, Hermiticity, six-face balance, charge, tail norm, escape
norm, bounded real-space-DG energy/density difference, and checkpoint
round-trip success.

**Step 4: Focused verification**

Run:

```bash
python3 tests/dg/run_dg_wannier_ground_state_mpi.py
python3 tests/dg/run_dg_buffer_wannier_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS.

**Step 5: Review, overlay build, and commit**

Resolve all Critical/Important findings and commit:

```bash
git commit -m "feat(dg): converge Wannier-coefficient ground state"
```

### Task 8: Propagate fixed Wannier orbitals through coefficient-only RT

**Prerequisite:** Tasks 5--7 passed.

**Files:**

- Create: `src/rt/dg/rt_dg_wannier_coefficients.f90`
- Modify: `src/rt/dg/CMakeLists.txt`
- Modify: the DG-only RT entry point selected by the existing route
- Create: `tests/dg/test_dg_wannier_coefficient_rt_mpi.f90`
- Create: `tests/dg/run_dg_wannier_coefficient_rt_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write RED coefficient-propagation tests**

Use a small generalized Hermitian system with known exponential solution.
Require:

- fixed Wannier orbitals;
- coefficient-only RK propagation;
- density reconstruction at every RK stage;
- Hartree/XC/vlocal update at every stage;
- external-field, nonlocal, and DG updates at every stage;
- overlap-metric norm conservation;
- restart equivalence;
- explicit failure on tail or window escape;
- no propagation of `nstate_frag` DC orbitals.

**Step 2: Implement the stage derivative**

Implement:

```fortran
subroutine evaluate_wannier_rt_stage(coefficients,time,dcoeff,diagnostics,ok,message)
  call reconstruct_stage_density(coefficients,...)
  call update_stage_potential(...)
  call update_stage_wannier_hamiltonian(...)
  call apply_overlap_inverse_hamiltonian(coefficients,dcoeff,...)
  dcoeff=(0d0,-1d0)*dcoeff
end subroutine
```

Every RK stage calls this routine with its stage coefficients.

**Step 3: Add fixed-space validity gates**

Monitor:

- \(C^\dagger S_WC\);
- particle number;
- energy behavior appropriate to the external field;
- tail norm;
- window escape norm;
- face balance;
- restart fingerprint.

Fail explicitly rather than updating the Wannier basis or falling back to
conventional RT.

**Step 4: Focused verification**

Run:

```bash
python3 tests/dg/run_dg_wannier_coefficient_rt_mpi.py
python3 tests/dg/run_dg_wannier_ground_state_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS.

**Step 5: Review, overlay build, and commit**

Resolve all Critical/Important findings, run a fresh overlay build including
the DG-only RT target, and commit:

```bash
git commit -m "feat(dg): propagate Wannier coefficients in real time"
```

## Final verification

After Task 8, and only if Task 5 passed, run:

```bash
python3 tests/dg/run_dg_buffer_window_projector.py
python3 tests/dg/run_dg_dc_buffer_core_faces_mpi.py
python3 tests/dg/run_dg_dc_direct_cg_mpi.py
python3 tests/dg/run_dg_dc_handoff_mpi.py
python3 tests/dg/run_dg_dc_direct_checkpoint_mpi.py
python3 tests/dg/run_dg_buffer_wannier_mpi.py
python3 tests/dg/run_dg_wannier_ground_state_mpi.py
python3 tests/dg/run_dg_wannier_coefficient_rt_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Then run:

- the normal DC Si64 reference through LCFO+EigenExa;
- the accepted projected real-space DG Si64 gate;
- Wannier GS restart equivalence;
- short and restart-split Wannier coefficient RT;
- a final fresh parent-prerequisite overlay build.

Use `@superpowers:finishing-a-development-branch` only after all required
verification passes.
