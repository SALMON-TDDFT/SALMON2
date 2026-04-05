# DG Hybrid RT Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a first hybrid MPI+OpenMP design for DG RT that covers both `H0` application and overlap solve.

**Architecture:** Keep MPI at fragment/group boundaries and add OpenMP only for rank-local DG RT kernels. Use state-batch partitioning as the first hybrid policy so output slices remain thread-safe without atomics. Keep BLAS/LAPACK threading disabled in the first implementation.

**Tech Stack:** Fortran, MPI, OpenMP, SALMON DG RT (`FF/FP/PP`, block operators, iterative overlap solve)

---

### Task 1: Confirm the active RT branch and preserve the current debug evidence

**Files:**
- Read: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90`
- Read: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Read: `/Users/otobetoshihito/Downloads/0/1/stdout.1.0`

**Step 1: Record the active derivative branch**

Confirm whether the large RT case enters:

- `apply_mixed_hamiltonian`
- `apply_matrix_blocks_batch`
- or dense `zgemm`

using the current derivative trace markers.

**Step 2: Save the result in a short debug note**

Write one short note under `docs/plans/` or append to the design doc with:

- active branch
- observed stop stage
- whether overlap solve was reached

**Step 3: Commit the note if it changes**

```bash
git add docs/plans/2026-03-26-dg-hybrid-rt-design.md
git commit -m "docs: record DG RT hybrid debug target"
```

### Task 2: Add a runtime switch for DG hybrid RT

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/io/salmon_global.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/io/inputoutput.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90`

**Step 1: Add input/runtime state**

Introduce one DG RT hybrid control, for example:

- `yn_dg_hybrid_rt`

with default disabled.

**Step 2: Broadcast and log it**

Make sure the flag is:

- read from namelist
- broadcast to ranks
- written to `variables.log`

**Step 3: Add per-fragment runtime state if needed**

If the apply/solve routines need a cached logical on `dg_frag`, add it to the DG fragment type.

**Step 4: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected:

- build succeeds

**Step 5: Commit**

```bash
git add src/io/salmon_global.f90 src/io/inputoutput.f90 src/rt/dg/rt_dg_fragment_types.f90
git commit -m "feat: add DG RT hybrid runtime switch"
```

### Task 3: Refactor `apply_matrix_blocks_batch` for explicit state-batch partitioning

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Test: small DG RT case under `/tmp`

**Step 1: Identify all writes to the output batch**

List the output slices updated inside `apply_matrix_blocks_batch`.

**Step 2: Introduce state-range helpers**

Add a small helper or inline logic to compute a thread-owned column range for `x(:, :)` and `y(:, :)`.

**Step 3: Add OpenMP over state ranges**

Use `!$omp parallel do` over thread-owned state batches so each thread updates only:

- `y(:, jbeg:jend)`

for its assigned states.

**Step 4: Keep block traversal unchanged**

Inside each thread, traverse the current block structure exactly as before so the math does not change.

**Step 5: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

**Step 6: Run a small RT smoke case**

Use a 1-rank and a small hybrid-friendly run if available.

Expected:

- no numerical regression
- no OpenMP race symptoms

**Step 7: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "feat: openmp-parallelize DG block apply by state batch"
```

### Task 4: Refactor `apply_mixed_hamiltonian` with the same state-batch policy

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Test: mixed DG+PW RT case under `/tmp`

**Step 1: Separate `FF`, `FP`, and `PP` contributions conceptually**

Document the contribution order inside the routine:

- fragment-fragment
- fragment-PW
- PW-fragment
- PW-PW

**Step 2: Apply the same state-batch threading**

Parallelize over state columns and ensure each thread owns a unique output slice.

**Step 3: Avoid nested threading**

If the routine calls other threaded kernels, keep only one OpenMP layer active.

**Step 4: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

**Step 5: Verify on a mixed DG+PW smoke case**

Expected:

- same output as serial or MPI-only baseline
- no hang inside `before-apply-mixed-h`

**Step 6: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "feat: openmp-parallelize DG mixed Hamiltonian apply"
```

### Task 5: Add OpenMP to `solve_overlap_operator_batch`

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Test: DG RT overlap solve cases under `/tmp`

**Step 1: Identify the local kernels inside the iteration**

Separate:

- overlap matvec
- local dot products
- local norms
- vector updates

**Step 2: Apply the same state-batch partition**

Use the same state-range ownership model as `H0` apply.

**Step 3: Use OpenMP reductions only for local scalar reductions**

Do not move MPI collectives inside threaded regions.

**Step 4: Keep iteration ordering unchanged**

The first hybrid version must not change the solver algorithm.

**Step 5: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

**Step 6: Verify on a small RT case**

Expected:

- same convergence behavior
- no NaN or deadlock

**Step 7: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "feat: openmp-parallelize DG overlap solve"
```

### Task 6: Add hybrid logging and execution guidance

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-26-dg-hybrid-rt-design.md`

**Step 1: Log hybrid status at runtime**

Add one concise runtime log line showing:

- hybrid enabled or disabled
- OpenMP thread count
- partition mode

**Step 2: Document thread policy**

Add a short note that BLAS threading should stay off for the first hybrid mode.

**Step 3: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_integrator_derivative.f90 src/rt/dg/rt_dg_fragment_ops.f90 docs/plans/2026-03-26-dg-hybrid-rt-design.md
git commit -m "docs: log DG RT hybrid runtime policy"
```

### Task 7: Verify hybrid behavior on representative runs

**Files:**
- Use existing runtime inputs and logs

**Step 1: Run serial or MPI-only baseline**

Capture reference output for a small DG RT case.

**Step 2: Run hybrid-enabled small case**

Use the same physics input with OpenMP threads enabled.

**Step 3: Compare outputs**

Check:

- coefficient norms
- observables
- absence of hangs at `H0`
- overlap solve completion

**Step 4: Run one larger hybrid smoke**

Use a reduced Fugaku-like layout if available.

**Step 5: Commit any test/log updates**

```bash
git add <relevant files>
git commit -m "test: verify DG RT hybrid H0 and overlap paths"
```
