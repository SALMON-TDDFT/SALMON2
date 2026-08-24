# WF+PW LCFO Divided-SCF Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace repeated full-cell Hybrid eigensolves by a DC-style fragment WF+PW SCF, followed by one distributed WF+PW LCFO diagonalization.

**Architecture:** Preserve the existing conventional DC+LCFO/Wannier90 route as the WF seed generator.  After the WF and windowed-PW complement are fixed, solve the Hybrid density fragment-by-fragment with the existing DC convergence controls, then assemble a distributed WF+PW LCFO generalized eigenproblem once.  Never gather full-system real-space orbitals or replicated full `H/S` matrices.

**Tech Stack:** Fortran 2008, SALMON DC/LCFO, MPI, BLAS/LAPACK, EigenExa or complex ScaLAPACK, CMake, Python source-contract tests, linked MPI fixtures.

---

Implementation must preserve all pre-existing user changes in this dirty worktree.  The
`tests/dg/run_*.py` and `tests/dg/check_*.py` programs in this worktree are
standalone tests and are not registered with CMake/CTest; execute them directly
and do not create a new test-registration layer.  Stage and commit only the
current task's new hunks.  For every file that was already modified before this
plan started, use `git add -p <file>` rather than `git add <file>`, then inspect
`git diff --cached` and reject every unrelated hunk before committing.  Use
`MPI=8` and `OMP_NUM_THREADS=1` for Si64; do not use a time-based cutoff.

### Task 1: Freeze the route contract

**Files:**
- Create: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`

**Step 1: Write the failing source-contract test**

Require one explicit divided route in `main_dft.f90` after the accepted W90/WPW basis exists.  Assert that the route calls `run_dg_hybrid_divided_scf`, then `assemble_dg_hybrid_lcfo_rows`, then exactly one generalized eigensolver call.  Assert that it does not call `run_dg_hybrid_self_consistent_ground_state` inside the divided branch and does not allocate an `ntarget`-by-`ntarget` replicated matrix.

```python
def test_divided_route_precedes_one_shot_lcfo():
    source = Path("src/gs/main_dft.f90").read_text()
    branch = source[source.index("if(yn_dg_hybrid_divided_scf=='y')"):]
    assert branch.index("run_dg_hybrid_divided_scf") < branch.index("assemble_dg_hybrid_lcfo_rows")
    assert branch.count("solve_dg_hybrid_generalized_") == 1
    assert "run_dg_hybrid_self_consistent_ground_state" not in branch.split("endif", 1)[0]
```

**Step 2: Run it and verify RED**

Run: `python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py`

Expected: FAIL because the divided branch and driver do not exist.

**Step 3: Keep the standalone RED test**

Keep the Python check beside the existing Hybrid source-contract tests.  Do not
register it with CMake and do not add production stubs merely to turn it green.

**Step 4: Commit the RED contract**

```bash
git add tests/dg/check_dg_hybrid_divided_lcfo_route.py
git commit -m "test(dg): specify divided WF+PW LCFO route"
```

### Task 2: Add fragment WF+PW catalog and ownership

**Files:**
- Create: `src/gs/dc/dg_hybrid_fragment_basis.f90`
- Create: `tests/dg/test_dg_hybrid_fragment_basis_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_fragment_basis_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing MPI test**

Build two fragments on two ranks with WF rows, projected PW rows, center-fragment IDs, and bounded buffer support.  Require deterministic fragment-local ordering (WF first, then PW), exactly-one center ownership, no off-fragment persistent values, and rejection of duplicate IDs or a support radius larger than the buffer.

**Step 2: Run it and verify RED**

Run: `python3 tests/dg/run_dg_hybrid_fragment_basis_mpi.py`

Expected: compile failure because `dg_hybrid_fragment_basis` is absent.

**Step 3: Implement the minimal catalog**

Expose a bounded type and constructor:

```fortran
type, public :: s_dg_hybrid_fragment_basis
  integer :: fragment_id=0, generation=0
  integer, allocatable :: global_ids(:), sector(:)
  complex(real64), allocatable :: buffer_values(:,:)
  integer(int64) :: provenance_fingerprint=0_int64
end type

subroutine build_dg_hybrid_fragment_basis(comm,fragment_id,wf_ids,wf_values,pw_ids,pw_values,&
    projector_radius,buffer_radius,basis,ok,message)
```

Validate collective metadata before allocation.  Reuse the accepted WF/PW projection and rank-filter outputs; do not recompute Wannier90 or invent a second PW selector.

**Step 4: Run focused tests GREEN**

Run: `python3 tests/dg/run_dg_hybrid_fragment_basis_mpi.py`

Expected: PASS for 1, 2, and 4 ranks with identical fingerprints.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_hybrid_fragment_basis.f90 tests/dg/test_dg_hybrid_fragment_basis_mpi.f90 tests/dg/run_dg_hybrid_fragment_basis_mpi.py
git add -p src/gs/dc/CMakeLists.txt
git diff --cached --check
git diff --cached
git commit -m "feat(dg): catalog fragment WF+PW bases"
```

### Task 3: Implement one divided Hybrid density step

**Files:**
- Create: `src/gs/dc/dg_hybrid_divided_scf.f90`
- Create: `tests/dg/test_dg_hybrid_divided_scf_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_divided_scf_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write a failing two-fragment fixture**

Use a small gapped block Hamiltonian with known full solution.  Test one callback-driven divided step: restrict the total potential, solve each local generalized problem, apply common occupations supplied by the existing DC callback, and return only unique-core density.  Check electron count and rank-decomposition invariance.  With an empty PW sector, require bitwise-identical control flow and tolerance inputs to the conventional DC callbacks.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py`

Expected: compile failure because the driver is absent.

**Step 3: Implement the callback interface**

```fortran
subroutine run_dg_hybrid_divided_scf(comm,global_point_count,core_ids,initial_density,&
    convergence_mode,threshold,update_total_potential,solve_fragments,assemble_core_density,&
    mix_dc_density,maximum_iterations,converged_density,iterations,convergence_value,ok,message)
```

The driver must use the same convergence quantity selected by SALMON (`rho_dne`, `norm_rho`, or `norm_rho_dng`) and call the existing DC mixer callback.  Do not add energy, boundary-density, or post-LCFO convergence gates.  Reject duplicate core ownership collectively.

**Step 4: Run GREEN and regress the existing Hybrid driver**

Run:

```bash
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/run_dg_hybrid_scf_mpi.py
```

Expected: PASS; the existing full-cell reference driver remains unchanged.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_hybrid_divided_scf.f90 tests/dg/test_dg_hybrid_divided_scf_mpi.f90 tests/dg/run_dg_hybrid_divided_scf_mpi.py
git add -p src/gs/dc/CMakeLists.txt
git diff --cached --check
git diff --cached
git commit -m "feat(dg): add divided Hybrid SCF driver"
```

### Task 4: Connect SALMON DC potential, occupations, mixing, and core density

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dcdft.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`
- Create: `tests/dg/check_dg_hybrid_divided_dc_controls.py`

**Step 1: Extend the failing contract**

Require a default-off `yn_dg_hybrid_divided_scf`, initialized and broadcast like the current Hybrid option.  Require the production adapter to pass the existing `convergence`, `threshold`, DC density mixer, and total `dc%rho_tot`; forbid new Hybrid density tolerances in this branch.

**Step 2: Run RED**

Run: `python3 tests/dg/check_dg_hybrid_divided_dc_controls.py`

Expected: FAIL because the input and adapter are absent.

**Step 3: Add the production adapter minimally**

Start the divided stage only after conventional DC convergence, accepted Wannier90 output, and fixed WPW construction.  Initialize it from `dc%rho_tot`.  Reuse DC potential update, occupation/chemical-potential, mixing, and core-density assembly routines through narrow callbacks or extracted shared helpers.  Do not duplicate their algebra in `main_dft.f90`.

**Step 4: Run contracts and existing route tests**

Run:

```bash
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: all PASS; old flags preserve old behavior.

**Step 5: Commit**

```bash
git add tests/dg/check_dg_hybrid_divided_lcfo_route.py tests/dg/check_dg_hybrid_divided_dc_controls.py
git add -p src/gs/main_dft.f90 src/gs/dc/dcdft.f90 src/io/salmon_global.f90 src/io/inputoutput.f90
git diff --cached --check
git diff --cached
git commit -m "feat(dg): connect divided Hybrid SCF to DC controls"
```

### Task 5: Assemble distributed WF+PW LCFO rows once

**Files:**
- Create: `src/gs/dc/dg_hybrid_lcfo.f90`
- Create: `tests/dg/test_dg_hybrid_lcfo_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_lcfo_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing MPI assembly test**

Construct a synthetic WF+PW basis whose supports cross one fragment face and include a nonlocal projector centered on the neighboring core.  Compare distributed row-owned `H/S` against a direct dense reference.  Verify exactly-once local, kinetic, nonlocal, and boundary contributions; confirm peak storage is row-local rather than replicated `M*M` per rank.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_lcfo_mpi.py`

Expected: compile failure because `assemble_dg_hybrid_lcfo_rows` is absent.

**Step 3: Implement bounded assembly**

```fortran
subroutine assemble_dg_hybrid_lcfo_rows(comm,bases,row_ids,apply_h,apply_s,&
    hrows,srows,peak_elements,operator_fingerprint,ok,message)
```

Reuse the cached bidirectional redistribution and current Hybrid `hpsi` callbacks.  Communicate only required neighboring support in the operator layer; the first correctness implementation may use the already accepted all-to-all map, but the persistent matrix representation must remain block/row distributed.  Never create `orbital_owned_full_values` or a replicated dense matrix.

**Step 4: Run GREEN and redistribution regressions**

Run:

```bash
python3 tests/dg/run_dg_hybrid_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_nonlocal_projector_range_mpi.py
```

Expected: PASS and decomposition-independent fingerprints.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_hybrid_lcfo.f90 tests/dg/test_dg_hybrid_lcfo_mpi.f90 tests/dg/run_dg_hybrid_lcfo_mpi.py
git add -p src/gs/dc/CMakeLists.txt
git diff --cached --check
git diff --cached
git commit -m "feat(dg): assemble distributed WF+PW LCFO rows"
```

### Task 6: Perform one final distributed eigensolve and publish the state

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_hybrid_generalized_eigensystem.f90`
- Modify: `src/gs/dc/dg_hybrid_ground_state_types.f90`
- Modify: `tests/dg/test_dg_hybrid_generalized_eigensystem_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`

**Step 1: Add failing one-shot assertions**

Count eigensolver invocations with a fixture callback and require exactly one after divided-SCF convergence.  Verify `HC-SCepsilon`, `C^HSC-I`, occupations, electron count, and collective failure without checkpoint publication.  Do not test post-LCFO density convergence because it is not part of DC+LCFO semantics.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py`

Expected: FAIL on the missing divided-LCFO publication path.

**Step 3: Wire the final solve**

Select the real EigenExa-compatible backend only for a validated real Gamma representation; otherwise retain the complex distributed ScaLAPACK solver.  Keep eigenvectors coefficient-row distributed and feed the existing occupied checkpoint writer in bounded tiles.  Remove no existing reference solver.

**Step 4: Run GREEN**

Run:

```bash
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_ground_state_types_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
```

Expected: all PASS and the invocation count is one.

**Step 5: Commit**

```bash
git add tests/dg/check_dg_hybrid_divided_lcfo_route.py
git add -p src/gs/main_dft.f90 src/gs/dc/dg_hybrid_generalized_eigensystem.f90 src/gs/dc/dg_hybrid_ground_state_types.f90 tests/dg/test_dg_hybrid_generalized_eigensystem_mpi.f90 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
git diff --cached --check
git diff --cached
git commit -m "feat(dg): finalize divided WF+PW LCFO state"
```

### Task 7: Si64 physical validation

**Files:**
- Create: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in`
- Create: `tests/dg/run_dg_hybrid_si64_divided_lcfo.py`
- Modify: `docs/plans/2026-08-24-wpw-lcfo-divided-scf-design.md`

**Step 1: Add the validation harness**

Run from an explicit input file, use eight MPI ranks and export `OMP_NUM_THREADS=1`.  Parse conventional DC convergence, Wannier90 acceptance, divided WF+PW convergence, exactly one final LCFO eigensolve, residuals, electron count, energy, gap, symmetry receipts, wall time, and peak memory.  Do not set `time_shutdown`, shell `timeout`, or an application time cutoff.

**Step 2: Run the focused standalone DG tests first**

Run:

```bash
python3 tests/dg/run_dg_hybrid_fragment_basis_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: all selected tests PASS.

**Step 3: Run Si64**

Run: `OMP_NUM_THREADS=1 python3 tests/dg/run_dg_hybrid_si64_divided_lcfo.py --mpi-ranks 8`

Expected: normal completion, DC density convergence under the existing input criterion, accepted W90/WPW basis, divided Hybrid convergence, one final LCFO diagonalization, valid checkpoint, and finite receipts.

**Step 4: Compare with the stored full-system Hybrid oracle**

Compare total energy, gap, occupied projector, density diagnostic, symmetry, memory, and timing.  Then repeat with the next supported buffer width and PW cutoff/packet count.  Record convergence; do not convert the diagnostic density difference into a production SCF gate.

**Step 5: Update the design evidence and commit**

```bash
git add tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in tests/dg/run_dg_hybrid_si64_divided_lcfo.py
git add -p docs/plans/2026-08-24-wpw-lcfo-divided-scf-design.md
git diff --cached --check
git diff --cached
git commit -m "test(dg): validate divided WF+PW LCFO on Si64"
```

### Task 8: Final review and regression gate

**Files:**
- Review all files changed by Tasks 1-7

**Step 1: Review ownership and scaling**

Confirm no full-system real-space orbital array, per-rank dense `H/S`, repeated final eigensolve, or SCF-time global occupied coefficient matrix was introduced.  Confirm buffer/projector support and neighboring-core contributions are complete.

**Step 2: Review physical semantics**

Confirm the conventional DC+LCFO/W90 seed route is unchanged, existing DC convergence settings are authoritative, metals remain rejected, and no post-LCFO density-SCF loop was added.

**Step 3: Run formatting and focused verification**

```bash
git diff --check
python3 tests/dg/run_dg_hybrid_fragment_basis_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_nonlocal_projector_range_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
OMP_NUM_THREADS=1 python3 tests/dg/run_dg_hybrid_si64_divided_lcfo.py --mpi-ranks 8
```

Expected: no diff errors, all focused tests PASS, and Si64 completes without a time cutoff.

**Step 4: Request code review**

Use `superpowers:requesting-code-review`.  Address findings one at a time under `superpowers:receiving-code-review`, rerunning the affected focused test after each change.

**Step 5: Commit review fixes only if needed**

```bash
git add <reviewed-files-only>
git commit -m "fix(dg): address divided WF+PW LCFO review"
```
