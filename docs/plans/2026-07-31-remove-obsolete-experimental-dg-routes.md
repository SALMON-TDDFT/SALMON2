# Obsolete Experimental DG Route Removal Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove every experimental DG execution path except the accepted buffer-supported overlapping-Wannier GS and V3-checkpoint-backed generalized-eigenvalue Exp coefficient RT, while preserving normal SALMON and normal DC LCFO plus EigenExa.

**Architecture:** Remove obsolete routes from their public input and dispatch boundary inward, using a source-contract inventory to prevent partial deletion. Move only genuinely shared retained dependencies into the overlapping-Wannier namespace; delete all other orphaned implementation, test, and sample code.

**Tech Stack:** Fortran 2008, MPI, CMake, Python source-contract and MPI fixture runners, ScaLAPACK, EigenExa.

---

### Task 1: Freeze the retained and forbidden route contract

**Files:**
- Create: `tests/dg/check_obsolete_dg_routes_removed.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Create: `docs/plans/2026-07-31-obsolete-dg-route-inventory.md`

**Step 1: Record the source inventory**

List every input variable, dispatch branch, CMake source, focused test, and
sample belonging to:

- direct real-space DG GS;
- DG-Fragment and Nodal RT;
- WPW GS and RT; and
- mixed-z, full-H seed, adaptive-DG-basis, and experimental DG promotion.

For every candidate, record its retained consumers from:

```bash
rg -n '<symbol>' src tests samples
cmake --build <existing-build> --target help
```

Mark normal LCFO, EigenExa, Wannier90, and overlapping-Wannier consumers as
protected.

**Step 2: Write the RED contract**

Create a Python checker with explicit forbidden tokens:

```python
FORBIDDEN_INPUTS = {
    "yn_dg_dc_local_periodic",
    "yn_dg_fragment_rt",
    "yn_dg_nodal_rt",
    "yn_dg_wpw_checkpoint_rt",
    "yn_dg_wpw_production",
}
FORBIDDEN_SOURCE_PREFIXES = (
    "src/common/dg_wpw_",
    "src/gs/dc/dg_wpw_",
    "src/rt/dg/rt_dg_wpw_",
    "src/rt/dg/rt_dg_nodal_",
)
REQUIRED_SOURCES = {
    "src/gs/dc/dg_overlapping_wannier_checkpoint.f90",
    "src/gs/dc/dg_overlapping_wannier_construction.f90",
    "src/rt/dg/rt_dg_overlapping_wannier.f90",
}
```

The checker must scan `src`, relevant `CMakeLists.txt`, `tests/dg`, and
`samples`; fail for every remaining forbidden entry; and fail if a required
retained source disappears.

Also assert that normal `lcfo.f90`, `eigen_subdiag_eigenexa.f90`, and their
CMake entries remain.

**Step 3: Run the RED contract**

Run:

```bash
python3 tests/dg/check_obsolete_dg_routes_removed.py
```

Expected: FAIL and print the complete current forbidden inventory.

**Step 4: Review the inventory**

Perform specification and code-quality reviews. Resolve every Critical and
Important classification error before deletion.

**Step 5: Commit**

```bash
git add tests/dg/check_obsolete_dg_routes_removed.py \
  tests/dg/check_dg_overlapping_wannier_route.py \
  docs/plans/2026-07-31-obsolete-dg-route-inventory.md
git commit -m "test(dg): inventory obsolete experimental routes"
```

### Task 2: Remove obsolete DG input and dispatch surfaces

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/rt/main_tddft.f90`
- Modify: `tests/dg/check_obsolete_dg_routes_removed.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Narrow the RED contract**

Add exact assertions that the four retained flags remain:

```python
REQUIRED_INPUTS = {
    "yn_dg_dc_overlapping_wannier",
    "yn_dg_overlapping_wannier_rt",
    "yn_dg_overlapping_wannier_rt_restart",
    "yn_dg_length_gauge",
}
```

Require every other inventory-classified DG route flag to be absent from
global declarations, namelists, defaults, broadcasts, logging, and input
validation.

Run the checker and confirm it still fails.

**Step 2: Remove obsolete global and namelist variables**

Delete the inventory-approved obsolete variables from:

- declarations in `salmon_global.f90`;
- `propagation`, `dc`, and `dg_fragment` namelists;
- default initialization;
- MPI broadcasts;
- `variables.log` output; and
- argument checks.

Do not change normal LCFO, EigenExa, Wannier90, or conventional RT options.

**Step 3: Remove obsolete dispatch**

Delete direct-DG and WPW branches from `main_dft.f90`. Retain only the
overlapping-Wannier branch and the unchanged normal branch.

Delete DG-Fragment, Nodal, and WPW dispatch and driver bodies from
`main_tddft.f90`. The coefficient-RT dispatch must remain before
conventional RT initialization and must return immediately.

**Step 4: Verify**

Run:

```bash
python3 tests/dg/check_obsolete_dg_routes_removed.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_dg_overlapping_wannier_metadata_mpi.py
python3 tests/dg/run_rt_dg_overlapping_wannier_mpi.py
git diff --check
```

Expected: route contract, metadata, and coefficient RT PASS; removal checker
still reports implementation files pending deletion.

**Step 5: Review and commit**

Review normal-route dispatch ordering, unknown-namelist failure behavior,
and forbidden fallback absence. Resolve all Critical/Important findings.

```bash
git add src/io/salmon_global.f90 src/io/inputoutput.f90 \
  src/gs/main_dft.f90 src/rt/main_tddft.f90 tests/dg
git commit -m "refactor(dg): remove obsolete route entry points"
```

### Task 3: Remove direct real-space DG ground-state implementation

**Files:**
- Delete: `src/common/dg_dc_direct_sipg.f90`
- Delete: `src/common/dg_ground_state_checkpoint.f90`
- Delete: `src/gs/dc/dg_dc_ground_state.f90`
- Delete: `src/gs/dc/dg_dc_ground_state_adapter.f90`
- Delete: `src/gs/dc/dg_dc_handoff.f90`
- Delete: `src/gs/dc/dg_dc_local_basis_ground_state.f90`
- Modify: `src/common/CMakeLists.txt`
- Modify: `src/gs/dc/CMakeLists.txt`
- Delete: direct-DG-only files listed in `docs/plans/2026-07-31-obsolete-dg-route-inventory.md`
- Modify: `tests/dg/check_obsolete_dg_routes_removed.py`

**Step 1: Add exact RED file assertions**

Make the removal checker fail while each inventory-approved direct-DG file
or CMake entry exists.

**Step 2: Prove retained code has no consumer**

For every file:

```bash
rg -n '<module-name>|<public-procedure>' src \
  --glob '!<file-being-removed>'
```

If overlapping-Wannier consumes a helper, move only that helper into the
smallest relevant `dg_overlapping_wannier_*.f90` module and add a focused
test before deletion.

**Step 3: Delete implementation and CMake entries**

Delete the approved files and remove their source-list entries. Do not
delete `dg_buffer_window_projector.f90`, `dg_generalized_algebra.f90`, or
`dg_dc_buffer_core_faces.f90` while retained consumers exist.

**Step 4: Verify**

Run the removal and retained route contracts, every overlapping-Wannier GS
fixture on 1/2/4/8 ranks, and `git diff --check`.

**Step 5: Review and commit**

Resolve all Critical/Important findings, then:

```bash
git add -A src/common src/gs/dc tests/dg
git commit -m "refactor(dg): remove direct real-space ground state"
```

### Task 4: Remove WPW ground-state and RT implementation

**Files:**
- Delete: `src/common/dg_wpw_*.f90`
- Delete: `src/gs/dc/dg_wpw_*.f90`
- Delete: `src/gs/dc/dg_wannier_pw_scf.f90`
- Delete: `src/rt/dg/rt_dg_wpw_*.f90`
- Modify: `src/common/CMakeLists.txt`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/rt/CMakeLists.txt`
- Delete: WPW-only tests and runners listed in the inventory
- Modify: `tests/dg/check_obsolete_dg_routes_removed.py`

**Step 1: Add RED assertions**

Require zero files matching the three WPW prefixes and zero WPW CMake
entries, flags, test runners, or samples.

**Step 2: Audit shared algebra**

Before deleting, use `rg` to prove whether
`dg_generalized_algebra.f90` is consumed by overlapping-Wannier. Retain it
only if it has a direct retained consumer; otherwise delete it in this task.

Do not preserve WPW compatibility wrappers.

**Step 3: Delete WPW code and focused assets**

Delete the implementation files, CMake entries, WPW-only fixtures,
checkers, documentation results, and samples recorded in the inventory.

**Step 4: Verify**

Run:

```bash
python3 tests/dg/check_obsolete_dg_routes_removed.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_dg_overlapping_wannier_operators_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_solver_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_scf_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_checkpoint_mpi.py
git diff --check
```

**Step 5: Review and commit**

Resolve all Critical/Important findings, then:

```bash
git add -A src tests samples docs
git commit -m "refactor(dg): remove obsolete WPW routes"
```

### Task 5: Remove DG-Fragment and Nodal RT implementation

**Files:**
- Delete: `src/rt/dg/rt_dg_fragment*.f90`
- Delete: `src/rt/dg/rt_dg_nodal_*.f90`
- Delete: `src/rt/dg/rt_dg_integrator_*.f90`
- Delete: legacy RT-DG support files identified as exclusive in the inventory
- Modify: `src/rt/CMakeLists.txt`
- Delete: DG-Fragment/Nodal-only tests and runners listed in the inventory
- Delete: `samples/exercise_dg_fragment_rt/`
- Modify: `tests/dg/check_obsolete_dg_routes_removed.py`

**Step 1: Add RED assertions**

Require the fragment, nodal, and obsolete integrator prefixes to be absent.
Explicitly protect:

```python
REQUIRED_SOURCES.add("src/rt/dg/rt_dg_overlapping_wannier.f90")
```

**Step 2: Classify residual RT-DG helpers**

For each remaining `src/rt/dg/rt_dg_*.f90`, record its consumers. Delete it
unless `rt_dg_overlapping_wannier.f90` directly requires it. Move any
small retained helper into `rt_dg_overlapping_wannier.f90`; do not retain a
generic legacy framework for one consumer.

**Step 3: Delete code, CMake entries, tests, and samples**

Remove all inventory-approved files and update CMake.

**Step 4: Verify**

Run the removal checker, retained route contract, coefficient RT fixture on
1/2/4/8 ranks, all overlapping-Wannier GS fixtures, and `git diff --check`.

**Step 5: Review and commit**

Review generalized-eigenvalue Exp isolation, restart identity, and absence
of conventional RT fallback. Resolve all Critical/Important findings.

```bash
git add -A src/rt tests/dg samples
git commit -m "refactor(dg): remove fragment and nodal RT"
```

### Task 6: Remove obsolete LCFO-to-DG export and experimental controls

**Files:**
- Modify: `src/gs/dc/lcfo.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Delete: obsolete export-only files and tests classified in the inventory
- Modify: `tests/dg/check_obsolete_dg_routes_removed.py`

**Step 1: Write RED assertions**

For every inventory-classified Flux/export symbol, assert absence only
after proving it has no normal LCFO, EigenExa, or Wannier90 consumer.

**Step 2: Remove DG-only export**

Delete DG seed/checkpoint export, WPW adapters, and experimental control
handling. Preserve normal LCFO matrix construction, diagonalization,
EigenExa dispatch, and normal outputs byte-for-byte where possible.

If `lcfo_flux.f90` has no retained consumer after prior tasks, delete it and
its CMake entry. Otherwise remove only DG-specific procedures.

**Step 3: Run normal regression first**

Run the unchanged normal-DC Si64 LCFO plus EigenExa fixture and compare:

- LCFO marker;
- EigenExa marker;
- converged energy;
- electron count; and
- expected output inventory.

Expected: unchanged within the existing regression tolerances.

**Step 4: Run retained DG verification**

Run every overlapping-Wannier fixture on 1/2/4/8 ranks and the removal
checker.

**Step 5: Review and commit**

Resolve all Critical/Important findings, then:

```bash
git add -A src/gs/dc src/io tests
git commit -m "refactor(dg): remove obsolete LCFO DG exports"
```

### Task 7: Remove stale tests, samples, and documentation

**Files:**
- Delete: obsolete DG files under `tests/dg/` listed in the inventory
- Delete: obsolete DG sample directories listed in the inventory
- Modify: `CMakeLists.txt`
- Modify: `tests/CMakeLists.txt`
- Modify: DG documentation that presents a removed route as executable
- Modify: `docs/plans/2026-07-31-obsolete-dg-route-inventory.md`
- Modify: `tests/dg/check_obsolete_dg_routes_removed.py`

**Step 1: Extend the RED checker**

Reject stale executable instructions, sample inputs using removed flags,
and test registrations naming removed executables. Historical plans may
describe removed work but must carry a clear historical/removed notice and
must not be linked as current instructions.

**Step 2: Delete obsolete assets**

Remove every approved file and update test registration. Do not delete
overlapping-Wannier tests, Si64 gate scripts, or normal regression assets.

**Step 3: Prove no orphan remains**

Run:

```bash
python3 tests/dg/check_obsolete_dg_routes_removed.py
rg -n 'yn_dg_dc_local_periodic|yn_dg_fragment_rt|yn_dg_nodal_rt|yn_dg_wpw_' \
  src tests samples CMakeLists.txt
```

Expected: checker PASS and `rg` finds no executable source/test/sample
reference.

**Step 4: Full focused verification**

Run all retained overlapping-Wannier runners on 1/2/4/8 ranks, normal
route contract tests, and `git diff --check`.

**Step 5: Review and commit**

Resolve all Critical/Important findings, then:

```bash
git add -A tests samples docs CMakeLists.txt
git commit -m "test(dg): remove obsolete route assets"
```

### Task 8: Production and clean-overlay acceptance

**Files:**
- Modify: `docs/plans/2026-07-31-obsolete-dg-route-inventory.md`
- Create: `docs/plans/2026-07-31-obsolete-dg-route-removal-results.md`

**Step 1: Run the final source contract**

Run:

```bash
python3 tests/dg/check_obsolete_dg_routes_removed.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
git diff --check
```

Expected: PASS.

**Step 2: Run the complete focused suite**

Run every retained overlapping-Wannier metadata, metric, construction,
operator, physical-matrix, solver, SCF, checkpoint, and coefficient-RT
runner on MPI 1, 2, 4, and 8 ranks.

**Step 3: Run clean-first parent-prerequisite overlay**

Build from `git archive HEAD` plus the current task diff and the documented
parent prerequisite overlay. Enable MPI, ScaLAPACK, and EigenExa. Run:

```bash
cmake --build <overlay-build> --clean-first -j1
```

Expected: SALMON links successfully with no removed source referenced.

**Step 4: Run normal and DG production gates**

Run:

1. unchanged normal-DC Si64 LCFO plus EigenExa;
2. fresh Si64 overlapping-Wannier GS checkpoint publication;
3. accepted checkpoint reuse with unchanged manifest hash;
4. short field-off and impulse coefficient RT; and
5. one-shot versus restart-split coefficient RT byte identity.

Reject logs containing direct-DG, WPW, DG-Fragment, Nodal, normal
checkpoint publication, or conventional RT markers from the retained
route.

**Step 5: Final reviews**

Perform specification and code-quality reviews over the complete branch.
Resolve every Critical and Important finding and repeat affected
verification.

**Step 6: Record and commit**

Record commands, hashes, runtimes, memory, normal regression comparison,
and review results.

```bash
git add docs/plans/2026-07-31-obsolete-dg-route-inventory.md \
  docs/plans/2026-07-31-obsolete-dg-route-removal-results.md
git commit -m "test(dg): verify obsolete route removal"
```

