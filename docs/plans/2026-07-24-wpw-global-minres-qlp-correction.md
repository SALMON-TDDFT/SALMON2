# WPW Global MINRES-QLP Correction Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a default-off fixed-H comparison route that solves the global projected correction equation with bounded-memory complex MINRES-QLP and tests it on the Task 16 B=6 gate.

**Architecture:** A new common solver reuses the existing compatible Ritz projectors, batched bounded H/S callbacks, global Gram callback, and operator snapshot validation. It applies a complex Hermitian Lanczos/MINRES-QLP recurrence independently to deterministic state batches, accepts only explicitly verified minimum-length S-orthogonal corrections, and plugs into the existing fixed-H callback selection without changing GMRES or normal production algebra.

**Tech Stack:** Fortran 2008, MPI collectives, LAPACK dense test oracle, existing SALMON bounded WPW callbacks, Python source-contract runners, CMake, EigenExa/ScaLAPACK overlay build.

---

## Execution boundary

- Work only in
  `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`
  on `codex/wpw-s-orthogonal-complement`.
- Expected starting HEAD is `5157d61` with clean status.
- Treat `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG` as read-only.
- Do not copy or commit its uncommitted prerequisites. Use a fresh `/tmp`
  parent-prerequisite overlay for full builds without `--delete`.
- Keep GMRES unchanged. Keep MINRES-QLP default `n`, mutually exclusive,
  fixed-H/continuation-only, and absent from normal SCF, checkpoint,
  publication, and RT.
- Use `@superpowers:test-driven-development` for Tasks 1 and 2,
  `@superpowers:systematic-debugging` for failures,
  `@superpowers:requesting-code-review` after every Task, and
  `@superpowers:verification-before-completion` before each commit.
- Resolve every Critical/Important review finding before proceeding.

### Task 1: Implement the complex MINRES-QLP kernel

**Files:**
- Create: `src/common/dg_wpw_global_minres_qlp_correction.f90`
- Modify: `src/common/CMakeLists.txt`
- Create: `tests/dg/test_dg_wpw_global_minres_qlp_correction_mpi.f90`
- Create: `tests/dg/run_dg_wpw_global_minres_qlp_correction_mpi.py`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`

**Step 1: Write the RED dense-oracle fixture**

Create a two-rank fixture with nonidentity SPD `S`, noncommuting Hermitian
`H`, S-orthonormal `Q`, and explicit

```fortran
P_R = I - Q*Q^H*S
P_L = I - S*Q*Q^H
A_i = P_L*(H-epsilon(i)*S)*P_R
b_i = -P_L*r_i
```

Cover nonsingular indefinite, compatible singular minimum-length,
incompatible singular, zero RHS, forced QLP transition, independent active
masks, iteration exhaustion, non-Hermitian callback corruption, stale
snapshot, nonfinite callback output, and scaled-equivalent cases. Use an
explicit S-complement basis plus `zgelss` pseudoinverse as the dense oracle.
Require:

```text
maxabs(Q^H*S*z) <= 1E-11
norm(b-A*z) <= max(tau_abs,tolerance*norm(b))
norm(z-z_pinv) <= scale-aware 1E-10
```

The incompatible singular case must fail with zero returned corrections.

**Step 2: Add the RED source/build contract**

In `tests/dg/test_dg_wpw_matrix_free_scf.py`, require the new source,
CMake entry, public solver/control/diagnostic types, Hermitian-probe
validation, QLP phase fields, explicit final residual, transactional
diagnostics, and absence of iteration-sized distributed storage.

Run:

```bash
python3 tests/dg/run_dg_wpw_global_minres_qlp_correction_mpi.py
python3 tests/dg/test_dg_wpw_matrix_free_scf.py
```

Expected: FAIL because the module and API do not exist.

**Step 3: Define the public API**

Implement:

```fortran
type,public :: s_dg_wpw_global_minres_qlp_controls
  integer :: max_iterations=64
  integer :: state_batch=8
  real(8) :: relative_tolerance=1d-2
  real(8) :: transition_condition=1d7
end type

type,public :: s_dg_wpw_global_minres_qlp_diagnostics
  real(8),allocatable :: initial_residual(:),recursive_residual(:)
  real(8),allocatable :: final_residual(:),s_orthogonality(:)
  real(8),allocatable :: equation_defect(:),projected_fraction(:)
  real(8),allocatable :: correction_norm(:),amplification(:)
  real(8),allocatable :: operator_norm(:),condition_estimate(:)
  real(8),allocatable :: hermitian_defect(:),lanczos_defect(:)
  integer,allocatable :: iterations(:),phase(:),transition_iteration(:)
  integer,allocatable :: numerical_rank(:),state_status(:)
  logical,allocatable :: compatible_singular(:),incompatible(:),converged(:)
end type
```

Expose solve and release procedures. Reuse the existing callback interfaces
or move only truly common interfaces/projector helpers into a small shared
module; do not change GMRES behavior.

**Step 4: Implement collective validation and Hermitian probes**

Collectively reject invalid controls, rank-disagreeing shapes/counts,
nonfinite inputs, malformed Q, allocation overflow, stale operator, and
rank-disagreeing masks/phases. Validate:

```text
max_iterations in [1,128]
state_batch in [1,16]
relative_tolerance finite in (0,1)
transition_condition finite and > 1
```

For deterministic finite probe columns, measure
`X^H*A*Y-(A*X)^H*Y` with a scale-aware threshold. Invoke the production
snapshot validator after every batched H/S action and after the final
explicit residual.

**Step 5: Implement bounded-memory complex MINRES-QLP**

Implement the Hermitian Lanczos three-term recurrence and stable complex
symmetric rotations from the approved design. Maintain only the fixed number
of distributed recurrence/solution vectors required by MINRES-QLP for each
state batch. Transition irreversibly to QLP when the condition estimate
exceeds the configured threshold or the triangular factor becomes
numerically rank deficient.

Zero inactive columns but preserve identical callback order. Reduce active
and phase masks with both `MPI_MIN` and `MPI_MAX` before continuing. Start
from zero, classify zero projected RHS without iteration, and compute the
minimum-length solution in compatible singular cases.

**Step 6: Implement explicit acceptance and transactional diagnostics**

Right-project the candidate, explicitly recompute `b-A*z`, and accept only:

```fortran
final_residual <= max(100d0*epsilon(1d0)*max(1d0,initial_residual), &
                      controls%relative_tolerance*initial_residual)
```

Also require finite output, S orthogonality, projection defect, consistent
snapshot, and no incompatible singular classification. A single failed state
zeros all output and returns failure after synchronizing diagnostics.

**Step 7: Run focused verification**

```bash
python3 tests/dg/run_dg_wpw_global_minres_qlp_correction_mpi.py
python3 tests/dg/run_dg_wpw_global_projected_correction_mpi.py
python3 tests/dg/test_dg_wpw_matrix_free_scf.py
python3 tests/dg/test_wpw_generalized_algebra.py
python3 tests/dg/run_dg_wpw_w_row_layout.py
git diff --check
```

Repeat the new MPI fixture with `-fcheck=all -fbacktrace`. Expected: PASS.

**Step 8: Review and commit**

Request a review focused on complex rotations, QLP recurrences,
minimum-length semantics, scale invariance, collective ordering, bounded
memory, and transactional failure. Resolve all Critical/Important findings,
rerun Step 7, then commit only Task 1 files:

```bash
git commit -m "feat(dg-wpw): add global MINRES-QLP correction solver"
```

### Task 2: Integrate the default-off fixed-H comparison route

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/common/dg_wpw_production_context.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

**Step 1: Add RED input and route contracts**

Require default, namelist, broadcast, variables log, finite/range validation,
fixed-H requirement, full mutual exclusion, route logging, fixed-H and every
continuation callback branch, and callback-free normal algebra.

Run:

```bash
python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py
python3 tests/dg/test_dg_wpw_matrix_free_scf.py
```

Expected: FAIL for missing MINRES-QLP controls and routing.

**Step 2: Add input controls**

Add:

```text
yn_dg_wpw_global_minres_qlp_correction = 'n'
dg_wpw_global_minres_qlp_max_iterations = 64
dg_wpw_global_minres_qlp_tolerance = 1D-2
dg_wpw_global_minres_qlp_state_batch = 8
dg_wpw_global_minres_qlp_transition_condition = 1D7
```

Validate values even while disabled. Reject simultaneous `y` with diagonal,
metric-block, fragment-local H-epsilon-S, global GMRES, or S-orthogonal
routes. Require fixed-H.

**Step 3: Store controls and diagnostics**

Use the production context's private allocatable polymorphic fields or add
equivalent private fields without making standalone context tests depend on
the new module. Allocate, overwrite safely, release on every cleanup path,
and never serialize them.

**Step 4: Implement the production callback**

Capture the current post-Ritz Q immediately before correction, snapshot the
completed bounded operator, call MINRES-QLP with the existing H/S/Gram and
validator callbacks, revalidate after return, and fail closed with zero
correction. Do not call GMRES or any local factor.

**Step 5: Route fixed-H and continuation**

Add a mutually exclusive callback branch in fixed-H, continuation trials, and
rollback. Log the effective route and controls once. At accepted application
counts 32/96/160 log all approved occupied/extra solver diagnostics. On
failure log available per-block maxima/counts before returning.

**Step 6: Extend runtime route tests**

Prove that the current Ritz observer is invoked immediately before the
MINRES-QLP callback, changes across outer updates, callback failure has no
fallback, normal algebra remains callback-free, and every continuation path
selects the same route.

**Step 7: Focused verification**

```bash
python3 tests/dg/run_dg_wpw_global_minres_qlp_correction_mpi.py
python3 tests/dg/run_dg_wpw_global_projected_correction_mpi.py
python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py
python3 tests/dg/test_dg_wpw_matrix_free_scf.py
python3 tests/dg/run_dg_wpw_production_operator_mpi.py
python3 tests/dg/test_wpw_generalized_algebra.py
python3 tests/dg/run_dg_wpw_occupied_w_basis_mpi.py
python3 tests/dg/run_dg_wpw_w_row_layout.py
git diff --check
```

Expected: PASS.

**Step 8: Review and commit**

Review route exclusivity, default-off scope, current-Q timing, continuation,
failure diagnostics, cleanup, and lack of serialization/publication changes.
Resolve all Critical/Important findings and rerun Step 7.

```bash
git commit -m "feat(dg-wpw): integrate global MINRES-QLP comparison route"
```

### Task 3: Parent-prerequisite overlay build and final review

**Files:**
- No intended source changes unless review finds a blocker.

**Step 1: Record boundaries**

Record branch, HEAD, clean status, recent commits, Task file scopes, merge
base, and a SHA-256 hash of the parent worktree porcelain status.

**Step 2: Create a fresh overlay**

Under a fresh `/tmp` directory, copy the parent excluding `.git`,
`.worktrees`, all `build*`, and `stage2d_wpw_runs`. Treat it read-only. Merge
the parent HEAD and this branch with `git merge-tree --write-tree`, apply the
parent's uncommitted patch and untracked prerequisites only to the temporary
source, and overlay only committed branch changes. Do not use `--delete`.

**Step 3: Configure and build**

```bash
cmake -S <overlay-source> -B <overlay-build> \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/mpifort \
  -DCMAKE_Fortran_FLAGS=-fallow-invalid-boz \
  -DUSE_MPI=ON -DUSE_SCALAPACK=ON -DUSE_EIGENEXA=ON -DUSE_WANNIER90=OFF
cmake --build <overlay-build> --target salmon -j1
```

Expected: `[100%] Built target salmon`.

**Step 4: Final review and verification**

Request cross-Task review. Resolve every Critical/Important finding. Rerun
all Task 2 focused tests, the overlay build, and `git diff --check`. Record
explicitly that the build used parent prerequisites absent from this branch.

### Task 4: Run and document the Task 16 B=6 gate

**Files:**
- Modify: `docs/plans/2026-07-24-wpw-global-minres-qlp-correction.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

**Step 1: Clone the restart**

Clone Task 16 B=6 into a unique fresh `/tmp` run. Preserve basis, normalized
seed, cutoff, tolerance, continuation, history, and publication settings.
Set every correction route `n` except global MINRES-QLP `y`; use 64,
`1E-2`, 8, and `1E7`.

**Step 2: Run**

Run the Task 3 overlay binary with `OMP_NUM_THREADS=1` and eight MPI ranks
through the fixed-H boundary.

**Step 3: Record evidence**

First require all 160 state corrections to pass and fixed-H to reach inner
32. If it does, record 32/96/160 occupied/extra outer residuals, solver
phase/transition/condition/rank/residual/norm diagnostics, search-metric
rank/discarded weights, Ritz defects, final info, and publication state.
Compare Task 16, Task 19, S-orthogonal, fragment H-epsilon-S, and global
GMRES. If it fails earlier, record the exact available failure diagnostics
and mark later fields as not reached.

**Step 4: Review and commit**

Review every transcription against raw logs. Resolve all Critical/Important
documentation findings and run `git diff --check`. Commit only the two docs:

```bash
git commit -m "docs(dg-wpw): record MINRES-QLP B=6 gate"
```

### Task 5: Decision checkpoint

Do not promote the route automatically.

Present the B=6 evidence for user review. Acceptance requires simultaneous
occupied/extra improvement with no solver nonconvergence, incompatible
singular state, search-rank collapse, Ritz inconsistency, or publication
regression. If rejected, keep both global routes default-off and use the
phase/rank/condition/state diagnostics to design an explicitly
state-partitioned policy.
