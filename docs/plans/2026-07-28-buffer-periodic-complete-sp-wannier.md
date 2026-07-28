# Buffer-Periodic Complete-s+p Wannier Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Derive the overlapping-Wannier target from the exact direct sum of the frozen occupied space and complete pseudopotential-backed atomic `s+p` residual shells, and pass the Si64 ground-state gate with the measured 384-function global basis.

**Architecture:** Reuse the tested SAWF pseudo-channel shell map and projection machinery without calling LCFO or Wannier90. Freeze the DC core-owned occupied projector, form its metric-orthogonal direct sum with every full-rank `s+p` projected residual direction in each periodic buffer box, localize only within that accepted space, and translate one representative basis to equivalent fragments. The manifest count controls shell completeness, while the measured direct-sum rank controls the retained basis count.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, SALMON DC pseudopotential/projector data, Python source-contract and Si64 evidence checkers.

---

### Task 1: Extract the complete-s+p manifest contract

**Files:**
- Create: `src/gs/dc/dg_overlapping_wannier_projection.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Test: `tests/dg/test_dg_overlapping_wannier_projection_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_projection_mpi.py`

**Step 1: Write the failing test**

Test atom-major channels for two atoms:

```fortran
call build_dg_complete_sp_manifest([1,2],channels,ok,message)
call require(ok.and.size(channels)==8,'two complete s+p shells')
call require(all(channels%l==[0,1,1,1,0,1,1,1]),'atom-major shell ordering')
```

Also reject incomplete counts and verify that the target is exactly four
times the core-owned atom count.

**Step 2: Verify RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_projection_mpi.py
```

Expected: compile failure because the module is absent.

**Step 3: Implement the manifest**

Move or expose only the generic channel type/count/map logic from
`lcfo_wannier_sawf.f90`.  Keep `lmax=1` explicit and do not introduce an
LCFO dependency.

**Step 4: Verify GREEN**

Run the focused test on 1, 2, 4, and 8 MPI ranks and run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
git diff --check
```

**Step 5: Review**

Perform specification and code-quality review.  Resolve all
Critical/Important findings before continuing.

### Task 2: Build pseudopotential-backed buffer projectors

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_projection.f90`
- Modify: `src/gs/main_dft.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_projection_mpi.f90`

**Step 1: Write the failing test**

Construct one `s` and three real `p` projectors in a periodic test box.
Require finite values, complete `(l,m)` coverage, periodic minimum-image
evaluation, and covariance under box translation.

**Step 2: Verify RED**

Run the focused projection fixture and confirm failure for the missing
projector builder.

**Step 3: Implement the builder**

Reuse the pseudo-channel real harmonics and the available pseudopotential
radial/projector tables.  Pass only the required immutable pseudopotential,
atom, lattice, grid, core, and buffer data from `main_dft.f90`.  Do not call
LCFO, EigenExa, WPW, or Wannier90.

**Step 4: Verify GREEN**

Run the projection fixture, route contract, and focused compilation.

**Step 5: Review**

Check periodic indexing, atom ownership, mixed-species ordering, temporary
allocation size, and pseudopotential fingerprint coverage.

### Task 3: Select the frozen-occupied plus complete-s+p target

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_projection.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_projection_mpi.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write the failing test**

Use a gauge-rotated candidate space with a two-dimensional occupied
projector and a complete four-channel manifest whose residual rank is
measured independently. Require:

- target rank equal to occupied rank plus complete-shell residual rank;
- exact occupied inclusion;
- invariant target projector under candidate gauge rotation;
- rejection when one `p` channel is linearly dependent;
- no arbitrary localizer-eigenvector padding.

**Step 2: Verify RED**

Run both focused fixtures and confirm failure at the new selection contract.

**Step 3: Implement selection**

Project the manifest into the candidate metric, remove the frozen occupied
component, diagonalize the projected complement Gram matrix, enforce a
scale-aware singular-value floor, and retain every independent residual
direction. Set the target to the occupied rank plus the measured residual
rank. Localize only by a unitary rotation inside this target.

**Step 4: Verify GREEN**

Run construction and projection fixtures on 1, 2, 4, and 8 ranks.

**Step 5: Review**

Check occupied inclusion, shell completeness, degeneracy handling,
candidate-gauge invariance, and bounded memory.

### Task 4: Connect core ownership and representative translation

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Test: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write the failing test**

Require manifest count `4 * core_owned_atom_count`, target count equal to
occupied rank plus measured complete-shell residual rank, correct fragment IDs
independent of MPI rank order, representative value/gradient translation,
and bond-center orbit closure after localization.

**Step 2: Verify RED**

Run the source contract and focused construction fixture.

**Step 3: Implement the production adapter**

Replace the numeric target-window authority with the measured direct-sum
rank. Retain the manifest count as the shell-completeness authority and the
candidate-window input only as the outer-space convergence control. Keep
complete buffer tails and the existing streaming closure.

**Step 4: Verify GREEN**

Run route, construction, metric, operator, and SCF focused checks.

**Step 5: Review**

Resolve all Critical/Important specification and quality findings.

### Task 5: Re-run the Si64 gate

**Files:**
- Modify: `tests/dg/data/si64_overlapping_wannier/inputfile_ow.in`
- Modify: `tests/dg/run_si64_overlapping_wannier_gate.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_gate.py`
- Modify: `docs/plans/2026-07-27-si64-overlapping-wannier-results.md`

**Step 1: Write the RED evidence requirements**

Require 32 complete-shell channels, measured residual rank 32, 48 retained
functions per fragment and 384 globally for 2x2x2 Si64,
complete-shell/projectability diagnostics, occupied rank 16 per fragment,
bond-center closure, and the existing full ground-state metrics.

**Step 2: Verify RED**

Run the checker on a fresh empty result root and confirm failure.

**Step 3: Run the minimum row**

Run 2x2x2, buffer 6, and the smaller candidate window.  Diagnose and fix
only root causes.  Do not increase to `s+p+d` and do not start Task 10.

**Step 4: Run the full matrix**

After the minimum row passes, run both decompositions, buffers 6 and 8, and
both candidate windows.  Record raw logs, memory, runtime, hashes, and route
checkpoint round trips.

**Step 5: Review and overlay build**

Perform specification review and code-quality review, resolve every
Critical/Important finding, and run a fresh clean-first
parent-prerequisite overlay build.

**Step 6: Commit**

Commit the accepted Task 9 implementation only after the Si64 gate passes:

```bash
git add src/gs/dc src/gs/main_dft.f90 src/io tests/dg docs/plans
git commit -m "test(dg): pass complete-sp Si64 Wannier gate"
```

Task 10 remains blocked unless this commit is backed by a passing Si64 gate.
