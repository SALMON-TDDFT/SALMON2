# DG Hybrid RT Design

**Goal:** Define a hybrid MPI+OpenMP design for the DG RT hot paths, with initial scope covering `H0` application and overlap solve on Fugaku-class runs.

## Context

Recent DG refactoring removed persistent mixed dense operators and moved normal RT paths to `FF/FP/PP` and block/apply forms. The current large-scale run now reaches RT step 1, but all ranks stop inside `calculate_time_derivative` before `after-h0`.

Observed runtime behavior from the current test logs:

- `momentum reduce` remains nontrivial during startup, with communication on the order of `6.7e-01 s`
- `hmat-reconstruct` at RT start is finite and no longer the dominant failure point
- all ranks reach `before-h0`
- no rank reaches `after-h0`

This means the current primary RT hot path is the `H0` application inside:

- `src/rt/dg/rt_dg_integrator_derivative.f90`
- `src/rt/dg/rt_dg_fragment_ops.f90`

and the next mandatory target is the overlap solve used in the same derivative path.

## Design Goals

1. Preserve the current MPI fragment/group decomposition.
2. Add OpenMP inside each rank for RT-local work.
3. Target both:
   - `H0` application
   - overlap solve
4. Avoid double-threading with BLAS/LAPACK in the first hybrid version.
5. Keep the design compatible with Fugaku execution where node-level hybrid tuning matters.

## Non-Goals

- No MPI topology redesign in the first step.
- No change to the mathematical algorithm of RT propagation.
- No immediate OpenMP expansion to every DG routine.
- No attempt to solve current RT hangs by changing physics or convergence criteria.

## Recommended Architecture

### MPI responsibilities

MPI remains responsible for:

- fragment/group decomposition
- neighbor communication
- halo/cache exchange
- collective metadata and block reductions

This keeps the existing weak-scaling-friendly decomposition intact.

### OpenMP responsibilities

OpenMP is introduced only for rank-local kernels:

- `apply_mixed_hamiltonian`
- `apply_matrix_blocks_batch`
- `solve_overlap_operator_batch`
- local vector operations inside the overlap iteration

The main parallel axis should be **state batch direction**, not block direction, for the first version. This minimizes write conflicts because each thread can own a disjoint `(:, istate_range)` slice of the output.

### BLAS threading policy

For the first hybrid version:

- keep BLAS/LAPACK internal threading disabled
- use OpenMP only in SALMON RT/DG kernels

This avoids oversubscription and unclear performance interactions.

## Candidate Approaches

### Approach 1: State-batch OpenMP

Parallelize over the occupied/state columns of `x(:, :)` and `y(:, :)`.

Pros:

- easiest to reason about
- naturally thread-safe for output slices
- applies to both `H0` apply and overlap solve

Cons:

- less ideal when `nstate_tot` is small
- may underuse threads if batch count is limited

### Approach 2: Block-parallel OpenMP with thread-local accumulation

Parallelize over matrix blocks and give each thread a private accumulation buffer before reduction.

Pros:

- matches block-structured DG operators
- potentially higher locality when block count is large

Cons:

- more memory
- more complex reductions
- riskier as a first hybrid step

### Approach 3: Rely on threaded BLAS

Leave Fortran loops mostly untouched and depend on threaded `matmul`/BLAS.

Pros:

- least code change

Cons:

- weak fit for current block/apply kernels
- poor control on Fugaku-style hybrid layouts
- higher risk of nested threading

## Recommendation

Use **Approach 1 first**, with explicit state-batch OpenMP in both `H0` apply and overlap solve.

This gives one coherent hybrid policy:

- MPI across fragments
- OpenMP across states within a rank
- BLAS threads off

It is the lowest-risk way to cover both target kernels from the start.

## Component Design

### 1. `H0` application

Target routines:

- `src/rt/dg/rt_dg_fragment_ops.f90`
  - `apply_mixed_hamiltonian`
  - `apply_matrix_blocks_batch`
- `src/rt/dg/rt_dg_integrator_derivative.f90`

Design:

- split state columns into thread-owned ranges
- each thread applies all required `FF/FP/PP` contributions for its state range
- no shared writes across threads for the same state range
- keep MPI/remote-coefficient fetch logic outside OpenMP regions or clearly single-threaded

Why this is preferred:

- avoids atomics on `y`
- avoids thread-private full-size `y_local`
- maps cleanly onto batched propagation

### 2. Overlap solve

Target routine:

- `src/rt/dg/rt_dg_fragment_ops.f90`
  - `solve_overlap_operator_batch`

Design:

- keep the Krylov/iterative structure unchanged
- OpenMP-parallelize:
  - overlap matvec
  - local dot products
  - local norms
  - vector updates
- keep MPI reduction only where global reductions are mathematically required

Important rule:

- matvec and local reductions should share the same state-batch partitioning policy as `H0` apply

This keeps the hybrid strategy uniform and easier to tune.

## Concurrency Rules

1. Do not call MPI collectives from inside multi-threaded regions in the first implementation.
2. Keep thread-private temporaries scoped per state batch.
3. Avoid OpenMP over block loops until state-batch mode is working and measured.
4. Treat remote coefficient pack/fetch phases as MPI-owned stages, not OpenMP-owned stages.

## Runtime Policy

The first implementation should be guarded by an explicit runtime switch, either:

- a new `&dg_fragment` input flag, or
- an environment variable such as `SALMON_DG_HYBRID_RT`

The runtime should log:

- whether DG hybrid RT is enabled
- OpenMP thread count
- whether BLAS threading is expected to be off
- chosen partition mode, initially `state-batch`

## Performance Expectations

Expected gains:

- better node-local utilization in `H0` apply
- better node-local utilization in overlap solve
- improved strong scaling at fixed node count by reducing per-rank serial work

Expected remaining limits:

- MPI communication and collectives remain
- startup momentum/block reductions remain separate bottlenecks
- if `nstate_tot` is too small, OpenMP efficiency may be limited

## Validation Plan

Validation should proceed in this order:

1. confirm the current `H0` branch actually taken in the large-scale run
2. implement state-batch OpenMP in the active `H0` kernel
3. verify numerical identity for small runs
4. apply the same pattern to overlap solve
5. measure hybrid scaling on representative RT cases

## Success Criteria

The design is successful if:

- RT step 1 progresses past the current `H0` bottleneck
- hybrid runs use OpenMP in `H0` apply and overlap solve without changing results
- node-level utilization improves with reduced MPI rank counts
- the implementation keeps the current `FF/FP/PP` and block-based structure intact
