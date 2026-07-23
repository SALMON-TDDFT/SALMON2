# WPW Global Projected Correction Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a default-off fixed-H comparison route that solves the global projected `H-epsilon S` correction equation with bounded-memory restarted GMRES and tests it on the Task 16 B=6 gate.

**Architecture:** A new common module applies compatible left/right Ritz projectors and a matrix-free projected correction operator through existing batched H/S and global Gram callbacks. It solves retained states in deterministic contiguous batches with restarted GMRES, returns only collectively validated S-orthogonal corrections, and exposes diagnostics to the existing LCFO fixed-H callback without changing normal production algebra.

**Tech Stack:** Fortran 2008, MPI collectives, existing SALMON bounded WPW H/S callbacks, Python source-contract runners, CMake, EigenExa/ScaLAPACK overlay build.

---

## Scope and invariants

- Work only in
  `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`
  on `codex/wpw-s-orthogonal-complement`.
- Do not copy, stage, or commit the parent worktree's uncommitted prerequisite implementation.
- Use a fresh `/tmp` parent-prerequisite overlay source for full builds, following the Task 2 procedure. Treat the parent worktree as read-only and do not use `--delete`.
- Keep the new route default `n`, fixed-H/continuation-only, and mutually exclusive with diagonal, metric-block, and fragment-local H-epsilon-S corrections.
- Do not change normal production algebra, checkpoint schema, publication, or RT handoff.
- Preserve the existing Ritz update, convergence, retained-search recurrence, and publication gates outside callback selection.
- Every Task ends with focused verification and code review. Resolve all Critical/Important findings before its commit.

### Task 1: Implement the projected operator and restarted GMRES kernel

**Files:**
- Create: `src/common/dg_wpw_global_projected_correction.f90`
- Modify: `src/common/CMakeLists.txt`
- Create: `tests/dg/test_dg_wpw_global_projected_correction_mpi.f90`
- Create: `tests/dg/run_dg_wpw_global_projected_correction_mpi.py`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`

**Step 1: Write the failing dense-oracle MPI fixture**

Create a two-rank fixture with a small distributed W/P layout, nonidentity SPD
`S`, noncommuting Hermitian `H`, two retained Ritz vectors, and three Ritz
values. Define test callbacks matching:

```fortran
subroutine apply_h_batch(context,xw,xp,yw,yp,info)
subroutine apply_s_batch(context,xw,xp,yw,yp,info)
subroutine global_gram_batch(left,right,nrow,nleft,nright,gram,info)
```

Require the public API:

```fortran
type(s_dg_wpw_global_correction_controls) :: controls
type(s_dg_wpw_global_correction_diagnostics) :: diagnostics

call solve_dg_wpw_global_projected_correction(context,MPI_COMM_WORLD,&
  apply_h,apply_s,global_gram,validate_operator_state,expected_epoch,&
  expected_fingerprint,qw,qp,eigenvalues,rw,rp,controls,zw,zp,&
  diagnostics,info)
```

The fixture must compare:

- `P_R x=x-Q(Q^dagger Sx)`;
- `P_L y=y-SQ(Q^dagger y)`;
- `P_L(H-epsilon*S)P_R x`;
- final correction against a dense KKT constrained solve, or equivalently a
  solve in an explicit S-complement basis; do not invert the singular
  full-space projected matrix;
- `maxabs(Q^dagger*S*z)<=1d-11`;
- explicit projected-equation residual against the design acceptance formula.

Add cases for zero projected RHS, happy breakdown, restart shorter than
`max_iter`, and a shortened final restart cycle.

**Step 2: Add RED source/build contracts**

In `tests/dg/test_dg_wpw_matrix_free_scf.py`, assert that the new source exists,
is listed in `src/common/CMakeLists.txt`, contains both projector formulas,
per-state diagnostics, explicit final residual recomputation, deterministic
state batching, and no global dense H/S assembly.

**Step 3: Run RED**

Run:

```bash
python3 tests/dg/run_dg_wpw_global_projected_correction_mpi.py
python3 tests/dg/test_dg_wpw_matrix_free_scf.py
```

Expected: FAIL because the module/API is absent.

**Step 4: Implement public types and collective validation**

Create:

```fortran
type,public :: s_dg_wpw_global_correction_controls
  integer :: restart=8
  integer :: max_iterations=32
  integer :: state_batch=8
  real(8) :: relative_tolerance=1d-2
end type

type,public :: s_dg_wpw_global_correction_diagnostics
  real(8),allocatable :: initial_residual(:),final_residual(:)
  real(8),allocatable :: relative_residual(:),s_orthogonality(:)
  real(8),allocatable :: equation_defect(:),projected_fraction(:)
  real(8),allocatable :: correction_norm(:),amplification(:)
  integer,allocatable :: iterations(:),restart_count(:)
  integer,allocatable :: breakdown_status(:),state_status(:)
  logical,allocatable :: zero_rhs(:),converged(:)
end type
```

Use `breakdown_status=0/1/2` for none/happy/failed breakdown and a named
`state_status` enumeration for zero RHS, converged, nonconverged, callback
failure, stale operator, and invalid result. Add
`release_dg_wpw_global_correction_diagnostics` and make repeated overwrite and
release safe.

Collectively reject restart outside `[1,16]`, maximum iterations outside
`[1,64]`, nonfinite tolerance outside `(0,1)`, state batch outside `[1,16]`,
rank-disagreeing state counts, malformed W/P/Q shapes, nonfinite inputs, and
overflow-prone allocation products.

**Step 5: Implement projectors and operator application**

Add private helpers equivalent to:

```fortran
right = x - Q * global_gram(Q, S*x)
left  = y - S*Q * global_gram(Q, y)
Ax    = P_L * (H - epsilon*S) * P_R * x
```

Cache `S*Q` only within one solve call. Reduce a collective error flag after
every H/S/Gram callback before any rank can return. The public solver receives
captured `expected_epoch`, `expected_fingerprint`, and a
`validate_operator_state(context,expected_epoch,expected_fingerprint,info)`
procedure. Invoke this validator collectively after every batched H and S
application, including projector applications and the final explicit residual
check. A stale result invalidates the solve before another Arnoldi update.

**Step 6: Implement bounded-memory restarted GMRES**

Process deterministic contiguous state batches of at most
`controls%state_batch`. Store at most `(restart+1)*state_batch` distributed
Krylov columns. Use modified Gram--Schmidt with one reorthogonalization pass,
complex Givens rotations, and per-state small Hessenberg systems.

Inactive or converged state columns must be zeroed while preserving identical
callback order on all ranks. Start every state from zero and every restart from
its explicitly recomputed residual.

After every provisional convergence or breakdown decision, reduce the integer
active mask with both `MPI_MIN` and `MPI_MAX`; reject disagreement before the
next callback. In the fixture, make `global_gram` deliberately return a
rank-local norm on a selected call while reporting `info=0`. This creates a
rank-local provisional convergence decision and verifies that the explicit
mask-agreement check fails collectively without deadlock.

Use:

```fortran
absolute_floor = 100d0*epsilon(1d0)*max(1d0,initial_norm)
target = max(absolute_floor,controls%relative_tolerance*initial_norm)
```

Accept zero RHS and happy breakdown only under the design criteria. Recompute
the final projected-equation residual explicitly; recursive estimates are
diagnostic only. Apply the right projector once more before returning and
validate `Q^dagger*S*z`.

Populate diagnostics transactionally in a candidate object throughout the
solve. On any validly initialized but failed solve, publish the candidate
diagnostics, including the failing state, breakdown/status, last explicit
residual, iteration count, correction norm, and amplification, while zeroing
all returned corrections. Only malformed inputs that prevent diagnostic
allocation may return an empty diagnostics object.

**Step 7: Add collective negative cases**

The MPI fixture must verify collective failure for rank-local shape mismatch,
nonfinite Ritz values/RHS/Q, invalid controls, callback failure, mid-cycle
epoch/fingerprint change, nonconvergence, and the injected rank-disagreeing
active mask described in Step 6. Confirm that all ranks return the same
nonzero `info` without deadlock and that initialized failures retain matching
diagnostics while corrections remain zero.

**Step 8: Run focused tests**

Run:

```bash
python3 tests/dg/run_dg_wpw_global_projected_correction_mpi.py
python3 tests/dg/test_dg_wpw_matrix_free_scf.py
python3 tests/dg/test_wpw_generalized_algebra.py
python3 tests/dg/run_dg_wpw_w_row_layout.py
git diff --check
```

Expected: PASS.

**Step 9: Review and commit**

Review projector signs, global reductions, GMRES least-squares updates,
breakdown handling, memory bounds, collective ordering, and explicit residual
acceptance. Resolve every Critical/Important finding and repeat Step 8.

Commit only the five Task 1 files:

```bash
git commit -m "feat(dg-wpw): add global projected correction solver"
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

Require:

- `yn_dg_wpw_global_projected_correction`, default `n`;
- namelist, broadcast, variables log, and y/n validation;
- controls with defaults restart `8`, max iterations `32`, tolerance `1d-2`,
  state batch `8`;
- exact range validation from the design;
- mutual exclusion across diagonal, metric-block, fragment-local H-epsilon-S,
  and global projected correction;
- fixed-H requirement;
- callback selection only in fixed-H and continuation;
- normal `wpw_algebra_step` remains callback-free;
- existing search-history propagation remains byte-for-byte present.

**Step 2: Run RED**

Run:

```bash
python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py
python3 tests/dg/test_dg_wpw_matrix_free_scf.py
```

Expected: FAIL for missing input and route plumbing.

**Step 3: Add input controls**

Declare, default, read, broadcast, log, and validate the new route and four
solver parameters. Keep the route default `n`. Reject invalid values even when
the route is disabled so variables logs and restarts cannot carry malformed
diagnostic state.

**Step 4: Store solver controls and diagnostics**

Store `s_dg_wpw_global_correction_controls` and the last diagnostics in the
production context. Do not add checkpoint serialization. Initialize controls
from input only when entering the fixed-H comparison boundary. Call
`release_dg_wpw_global_correction_diagnostics` next to
`release_dg_wpw_h_epsilon_s_factor` in the existing cleanup block in
`src/gs/dc/lcfo_flux.f90`. Test repeated callback overwrite, failed-solve
diagnostic publication, and final cleanup.

**Step 5: Add the LCFO callback**

Add `wpw_global_projected_precondition` with the existing
`dg_wpw_preconditioner` signature. It must:

1. snapshot the bounded operator epoch and layout fingerprint;
2. invoke `solve_dg_wpw_global_projected_correction` with current `Q`, Ritz
   values, residuals, completed fixed-H H/S callbacks, global Gram, captured
   identifiers, and a production validator callback;
3. let the solver collectively revalidate epoch/fingerprint after every H/S
   batch and once more after return;
4. return the globally S-orthogonal correction without extra damping;
5. fail closed with no fallback.

Do not call the fragment-local H-epsilon-S factor from this route.

**Step 6: Route fixed-H and continuation**

Add one mutually exclusive callback branch after metric-block and
fragment-local H-epsilon-S selection and before diagonal selection. Use the
same branch in every continuation trial and rollback. Since no factor is
cached, do not add lambda-refresh logic or preserve Krylov vectors.

Log the effective route and solver controls once. At callback counts
32/96/160, log the occupied/extra maxima and counts required by the design.
Keep existing window-state, search-metric, and Ritz diagnostics unchanged.

**Step 7: Add lifecycle and default-off runtime coverage**

Extend the two-rank matrix-free fixture or a dedicated LCFO source contract to
verify:

- default controls leave the existing no-callback route unchanged;
- explicit global route invokes the new callback;
- stale operator epoch/layout is rejected collectively;
- continuation starts a new zero-initialized solve after lambda changes;
- failed GMRES does not publish a correction or fall back.

**Step 8: Run focused and regression tests**

Run:

```bash
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

**Step 9: Review and commit**

Review default-off compatibility, input validation, mutual exclusion,
fixed-H-only routing, callback collective order, lifecycle, diagnostics,
normal-route invariance, and cleanup. Resolve every Critical/Important finding
and repeat Step 8.

Commit only the seven Task 2 files:

```bash
git commit -m "feat(dg-wpw): integrate global projected correction route"
```

### Task 3: Full parent-prerequisite overlay build and final code review

**Files:**
- No intended source changes unless review finds a Critical/Important issue.

**Step 1: Verify the worktree boundary**

Record branch, HEAD, clean status, recent commits, and a hash of the parent
worktree's porcelain status. Confirm Tasks 1 and 2 committed only their listed
files.

**Step 2: Create a fresh temporary source**

Under `/tmp`, copy the parent worktree while excluding `.git`, `.worktrees`,
all `build*` directories, and
`samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_runs`.
Treat the parent as read-only.

Apply `merge-base..HEAD` as diffs rather than blindly replacing files that
also contain parent prerequisite edits. Overlay only the committed branch
changes. Do not use `--delete`.

**Step 3: Configure and build**

Configure a new build with:

```bash
cmake -S <overlay-source> -B <overlay-build> \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/mpifort \
  -DCMAKE_Fortran_FLAGS=-fallow-invalid-boz \
  -DUSE_MPI=ON -DUSE_SCALAPACK=ON -DUSE_EIGENEXA=ON -DUSE_WANNIER90=OFF
cmake --build <overlay-build> --target salmon -j1
```

Expected: `[100%] Built target salmon`.

**Step 4: Run final review and verification**

Request review across Tasks 1 and 2. Resolve every Critical/Important finding.
Then rerun all Task 2 focused tests, the full overlay build, and
`git diff --check`.

If fixes are required, amend only the appropriate Task commit or create a
small scoped fix commit. Record explicitly that the build depended on
uncommitted parent prerequisites not included in this branch.

### Task 4: Run the Task 16 B=6 physical gate

**Files:**
- Modify: `docs/plans/2026-07-23-wpw-global-projected-correction.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

**Step 1: Clone the comparison input**

Clone the Task 16 B=6 restart into a fresh `/tmp` run directory. Preserve its
basis, normalized seed, cutoff, tolerance, continuation, history, and
publication settings.

Set:

```text
yn_dg_wpw_preconditioner = 'n'
yn_dg_wpw_metric_preconditioner = 'n'
yn_dg_wpw_h_epsilon_s_correction = 'n'
yn_dg_wpw_s_orthogonal_pw = 'n'
yn_dg_wpw_global_projected_correction = 'y'
```

Keep comparison controls restart `8`, max iterations `32`, tolerance `1d-2`,
and state batch `8`.

**Step 2: Run eight MPI ranks**

Run the Task 3 overlay binary through the fixed-H boundary with
`OMP_NUM_THREADS=1` and eight MPI ranks.

**Step 3: Record evidence**

At inner 32/96/160 record all design diagnostics, outer occupied/extra
residuals, search-metric effective rank/discarded weights, Ritz defects, final
`info`, and publication state.

Compare against:

- Task 16 no-precondition baseline;
- Task 19 metric-block correction;
- S-orthogonal complement gate;
- rejected fragment-local H-epsilon-S gate.

Interpret improvement only if both occupied and extra residuals improve
without GMRES nonconvergence, rank collapse, Ritz inconsistency, or
publication regression.

**Step 4: Review, verify, and commit results**

Review all transcribed values and ratios against the raw logs. Run
`git diff --check`. Resolve every Critical/Important documentation finding.

Commit only the two result documents:

```bash
git commit -m "docs(dg-wpw): record global correction B=6 gate"
```

#### Task 4 execution result

The Task 16 restart was cloned to
`/tmp/20260723_task4_global_projected_b6` and run with the Task 3 overlay
binary on eight MPI ranks with `OMP_NUM_THREADS=1`. The copied input retained
the Task 16 basis, normalized seed, cutoff, tolerance, continuation, search
history, and publication settings. The variables log confirms that diagonal,
metric-block, fragment-local H-epsilon-S, and S-orthogonal-complement routes
were all `n`; only global projected correction was `y`, with restart `8`,
maximum iterations `32`, tolerance `1.0E-2`, and state batch `8`.

The physical prerequisites reproduced: DC-SCF converged in 87 iterations,
the occupied-W tail ratio was `7.2686E-03` (warning), and the density seed had
captured norm `1`, projection residual `7.7960E-11`, normalization residual
`1.3946E-15`, projected charge error `5.9380E-11`, and source-to-DC residual
`6.0445E-02` (warning).

The first fixed-H correction failed closed before completing inner iteration
1. Seven of 160 state solves did not reach the requested relative tolerance
within 32 GMRES iterations. The failure summary reported
`info=1`, `post_snapshot_info=0`, maximum final projected-equation residual
and equation defect both `4.6571063211358618E-04`, maximum iterations `32`,
`nonconverged=7`, and `failed_breakdown=0`. Because no correction application
was accepted, the inner 32/96/160 design diagnostics, outer occupied/extra
residuals, search-metric rank/discarded weights, and Ritz post/direct defects
do not exist for this run; reporting baseline values in those slots would
misrepresent them as global-route measurements. The run instead emitted
`[DG-WPW-LOCAL-FAIL] fixed_h_stage=algebra iter=1 info=1`, performed the
existing full-LCFO fallback, exited normally, and published no WPW
checkpoint/manifest/RT state.

This gate is rejected. Unlike Task 16 no-precondition, Task 19 metric-block,
the S-orthogonal complement, and the rejected fragment-local H-epsilon-S
route, the global route did not reach inner 32 and therefore cannot show
simultaneous occupied/extra improvement. It also violates the explicit
no-GMRES-nonconvergence acceptance condition. The route remains default-off
and fixed-H/continuation-only; it is not promoted to normal SCF, checkpoint,
publication, or RT.

### Task 5: Decision checkpoint

Do not promote this route into normal outer SCF, checkpoint schema,
publication, or RT handoff automatically.

Present the B=6 evidence for user review. If accepted, create a separate
promotion design with restart/provenance semantics. If rejected, keep the
route default-off and use GMRES convergence/residual evidence to decide
between MINRES-QLP and an explicitly state-partitioned policy.
