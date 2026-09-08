# Si64 Hybrid Self-Consistent Ground-State Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build and validate a self-consistent Si64 ground state in the fixed Wannier plus windowed-PW hybrid basis, first with complex ScaLAPACK and then with an adaptively stopped block-CG solver.

**Architecture:** Keep the hybrid basis and metric immutable during the previously validated bounded two-point Anderson density loop (simple mixing plus one residual-history correction, not full Pulay).  Rebuild only the density-dependent sparse Hamiltonian, solve its occupied generalized eigenspace, reconstruct the full-cell density, and publish a multi-orbital checkpoint only after collective physical gates pass.  Use a complex Cholesky plus ScaLAPACK `PZHEEVD` path as the Si64 reference before enabling the matrix-free block-CG replacement.

**Tech Stack:** Fortran 2008, MPI, ScaLAPACK/BLACS, LAPACK test oracles, SALMON density/potential and mixing modules, CMake, Python MPI runners.

---

### Task 1: Freeze the revised route contract

**Files:**
- Modify: `docs/plans/2026-08-20-windowed-pw-wannier-hybrid.md`
- Create: `tests/dg/check_dg_hybrid_self_consistent_route.py`

**Step 1: Write the failing route test**

Require this order in the eventual production adapter:

```text
fixed hybrid basis -> initial DC+LCFO density -> update potential
-> sparse H rebuild -> occupied generalized solve -> density reconstruction
-> bounded two-point Anderson update -> convergence gates -> multi-orbital checkpoint
```

Reject a one-shot checkpoint publication, basis reselection inside SCF, the
single-vector RT checkpoint on this route, or RT dispatch before SCF convergence.

**Step 2: Run the route test and verify RED**

Run: `python3 tests/dg/check_dg_hybrid_self_consistent_route.py`

Expected: FAIL because no hybrid SCF production adapter exists.

**Step 3: Amend the old Task 10/11 wording**

Mark the Si8 one-shot generalized-spectrum route as superseded by the approved
Si64 SCF design.  Keep Si8 only as an algebraic fixture.

**Step 4: Commit**

```bash
git add docs/plans/2026-08-20-windowed-pw-wannier-hybrid.md \
  tests/dg/check_dg_hybrid_self_consistent_route.py
git commit -m "test: require self-consistent hybrid ground state"
```

### Task 2: Add a distributed multi-orbital state contract

**Files:**
- Create: `src/gs/dc/dg_hybrid_ground_state_types.f90`
- Create: `tests/dg/test_dg_hybrid_ground_state_types_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_ground_state_types_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing MPI tests**

Cover row-owned `coefficients(nowned,noccupied)`, occupations, eigenvalues,
basis/metric/operator provenance, exactly-once row ownership, finite payloads,
electron count, and decomposition-independent fingerprints.  Add REDs for a
missing row, duplicate row on one rank, negative occupation, wrong electron
count, stale metric receipt, and rank-disagreeing dimensions.

**Step 2: Run MPI 1/2/4/8 and verify RED**

Run: `python3 tests/dg/run_dg_hybrid_ground_state_types_mpi.py`

Expected: compile failure because the module is absent.

**Step 3: Implement the minimal contract**

Define a ground-state result type with allocatable distributed coefficient
matrix, occupations, eigenvalues, convergence receipts, and a complete provenance
fingerprint.  Use checked int64 extent/byte arithmetic and collective allocation
cleanup.

**Step 4: Run MPI 1/2/4/8 and commit**

Expected: PASS with identical fingerprints.

```bash
git add src/gs/dc/dg_hybrid_ground_state_types.f90 \
  tests/dg/test_dg_hybrid_ground_state_types_mpi.f90 \
  tests/dg/run_dg_hybrid_ground_state_types_mpi.py CMakeLists.txt
git commit -m "feat: define distributed hybrid occupied state"
```

### Task 3: Implement complex generalized eigensystem validation

**Files:**
- Create: `src/gs/dc/dg_hybrid_generalized_eigensystem.f90`
- Create: `tests/dg/test_dg_hybrid_generalized_eigensystem_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write the dense-oracle tests**

Use complex Hermitian `H`, positive-semidefinite rank-revealed `S`, a degenerate
occupied block, and fractional occupations.  Verify `H C-S C epsilon`,
`C^H S C-I`, occupied projector invariance under degenerate rotations, electron
count, and MPI decomposition independence.  Add singular-active-space,
non-Hermitian, nonfinite, and stale-provenance REDs.

**Step 2: Verify RED**

Run: `python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py`

**Step 3: Implement the complex ScaLAPACK reference adapter**

Use distributed complex Cholesky reduction, `PZHEEVD`, and back transformation.
Do not route the complex hybrid pencil through the real-only EigenExa adapter and
do not replicate full matrices.  Return occupied projector and residual receipts,
not a gauge-dependent eigenvector comparison.  Focused tests may use LAPACK
`ZHEGV` only as an independent oracle.

**Step 4: Run MPI 1/2/4/8 and commit**

```bash
git add src/gs/dc/dg_hybrid_generalized_eigensystem.f90 \
  tests/dg/test_dg_hybrid_generalized_eigensystem_mpi.f90 \
  tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py CMakeLists.txt
git commit -m "feat: validate hybrid generalized occupied states"
```

### Task 4: Reconstruct density from distributed occupied states

**Files:**
- Create: `src/gs/dc/dg_hybrid_density.f90`
- Create: `tests/dg/test_dg_hybrid_density_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_density_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing density tests**

Materialize bounded basis tiles and reconstruct
`rho(r)=sum_a occupation(a)*|sum_i basis_i(r) C_ia|^2` without retaining the
full basis or all spatial orbitals.  Compare with a dense complex oracle.  Test
degenerate occupied rotations, fractional occupations, zero-local-point ranks,
electron-number quadrature, nonfinite callbacks, and allocation cleanup.

**Step 2: Verify RED**

Run: `python3 tests/dg/run_dg_hybrid_density_mpi.py`

**Step 3: Implement bounded tile accumulation**

Batch occupied columns and spatial points.  Reuse the immutable basis-provider
receipt, use checked extents, and expose density and particle-number fingerprints.

**Step 4: Run MPI 1/2/4/8 and commit**

```bash
git add src/gs/dc/dg_hybrid_density.f90 tests/dg/test_dg_hybrid_density_mpi.f90 \
  tests/dg/run_dg_hybrid_density_mpi.py CMakeLists.txt
git commit -m "feat: reconstruct hybrid occupied density"
```

### Task 5: Add the nonlinear SCF controller with validated density-history mixing

**Files:**
- Create: `src/gs/dc/dg_hybrid_scf.f90`
- Create: `tests/dg/test_dg_hybrid_scf_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_scf_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write a synthetic nonlinear RED/GREEN fixture**

Use a small density-dependent Hermitian Hamiltonian with a known fixed point.
Require convergence from a DC-like initial density using the previously validated
two-point Anderson state.  Add fixtures for simple-mixing startup, bounded acceleration,
oscillation detection, history rejection/reset, decreasing mixing factor,
inner-solve failure, electron-count failure, and no checkpoint on failure.

**Step 2: Verify RED**

Run: `python3 tests/dg/run_dg_hybrid_scf_mpi.py`

**Step 3: Implement callback-driven outer SCF**

The controller owns convergence policy but calls the existing bounded two-point
Anderson mixer.  Callbacks update the potential, assemble sparse `H`, solve occupied
states, and reconstruct output density.  Keep the basis and `S` fingerprints
constant throughout the loop.

**Step 4: Run MPI 1/2/4/8 and commit**

```bash
git add src/gs/dc/dg_hybrid_scf.f90 tests/dg/test_dg_hybrid_scf_mpi.f90 \
  tests/dg/run_dg_hybrid_scf_mpi.py CMakeLists.txt
git commit -m "feat: converge hybrid density with Pulay mixing"
```

### Task 6: Implement adaptive block-CG stopping independently

**Files:**
- Create: `src/gs/dc/dg_hybrid_block_cg.f90`
- Create: `tests/dg/test_dg_hybrid_block_cg_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_block_cg_mpi.py`
- Modify: `CMakeLists.txt`

**Step 1: Write failing adaptive-control tests**

Test a convergent generalized problem, plateau, monotone growth, alternating
growth, clustered/degenerate occupied states, rank-deficient metric, and warm
start.  Assert a low initial cap, adaptive target derived from outer density
residual, explicit stop reason, exact iteration count, and no unconditional
100-iteration behavior.

**Step 2: Verify RED**

Run: `python3 tests/dg/run_dg_hybrid_block_cg_mpi.py`

**Step 3: Implement matrix-free block-CG**

Reuse the sparse exchange plans and metric solver.  Reorthogonalize in `S`, lock
only converged invariant subspaces, and compare occupied projectors across warm
starts.  Do not admit a failed inner solve to the outer SCF.

**Step 4: Run MPI 1/2/4/8 and commit**

```bash
git add src/gs/dc/dg_hybrid_block_cg.f90 tests/dg/test_dg_hybrid_block_cg_mpi.f90 \
  tests/dg/run_dg_hybrid_block_cg_mpi.py CMakeLists.txt
git commit -m "feat: solve hybrid occupied space with adaptive block CG"
```

### Task 7: Compare ScaLAPACK and block-CG inside the same SCF

**Files:**
- Modify: `src/gs/dc/dg_hybrid_scf.f90`
- Modify: `tests/dg/test_dg_hybrid_scf_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_scf_mpi.py`

**Step 1: Write a failing paired-solver test**

Run the same nonlinear fixture with complex ScaLAPACK and adaptive block-CG.  Compare the
converged occupied projector, density, energy, electron count, and symmetry
receipt.  Require fewer inner iterations after warm starts and verify that
intentional over-solving is rejected by policy.

**Step 2: Verify RED, implement solver selection, verify GREEN**

Run MPI 1/2/4/8.  Expected final projector/density agreement within declared
tolerances and identical physical fingerprints.

**Step 3: Commit**

```bash
git add src/gs/dc/dg_hybrid_scf.f90 tests/dg/test_dg_hybrid_scf_mpi.f90 \
  tests/dg/run_dg_hybrid_scf_mpi.py
git commit -m "test: compare hybrid EigenExa and block CG SCF"
```

### Task 8: Extend checkpointing to the occupied manifold

**Files:**
- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`
- Modify: `tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`

**Step 1: Write failing multi-orbital round-trip tests**

Write with one MPI decomposition and read with 1/2/4/8 ranks.  Verify all
occupied coefficients, occupations, eigenvalues, SCF convergence receipts, and
complete basis/operator provenance.  Add stale-SCF, incomplete-state, changed
occupation, and old single-vector-format REDs.

**Step 2: Verify RED**

Run: `python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`

**Step 3: Version and implement the new format**

Keep atomic publication and hostile-extent checks.  Do not silently reinterpret
the existing single-vector version.

**Step 4: Run MPI 1/2/4/8 and commit**

```bash
git add src/rt/dg/rt_dg_hybrid_checkpoint.f90 \
  tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90 \
  tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
git commit -m "feat: checkpoint hybrid occupied manifold"
```

### Task 9: Connect the diagnostic Si64 complex ScaLAPACK SCF route

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/check_dg_hybrid_self_consistent_route.py`
- Create: `tests/dg/run_dg_hybrid_si64_scf.py`
- Create: `docs/verification/si64-hybrid-self-consistent-ground-state.md`

**Step 1: Extend the route RED**

Require an explicit default-off diagnostic flag, fixed basis receipts across
iterations, existing Pulay initialization, complex ScaLAPACK reference selection, and
checkpoint publication only after all convergence gates.

**Step 2: Implement the minimal production adapter**

Use the converged DC+LCFO density as iteration zero.  Reuse SALMON potential
updates and mixing.  Require a ScaLAPACK-enabled build; do not silently fall back
to real EigenExa.  Do not enter RT or block-CG in this task.

**Step 3: Run focused tests before the material job**

```bash
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
python3 tests/dg/run_dg_hybrid_scf_mpi.py
python3 tests/dg/run_dg_hybrid_density_mpi.py
git diff --check
```

Expected: PASS.

**Step 4: Run Si64 with agreed resources**

Use the existing Si64 input, unchanged MPI rank count, and
`OMP_NUM_THREADS=1`.  Record every outer density/energy residual, ScaLAPACK
residual, Pulay action, electron count, symmetry defect, elapsed time, and
per-rank peak RSS.  Do not start RT.

**Step 5: Review results before enabling block-CG**

The review gate rejects unexplained oscillation, basis provenance changes,
nonmonotone long-period cycling, or checkpoint publication before convergence.

**Step 6: Commit**

```bash
git add src/gs/main_dft.f90 src/io/salmon_global.f90 src/io/inputoutput.f90 \
  tests/dg/check_dg_hybrid_self_consistent_route.py \
  tests/dg/run_dg_hybrid_si64_scf.py \
  docs/verification/si64-hybrid-self-consistent-ground-state.md
git commit -m "feat: converge Si64 hybrid reference ground state"
```

### Task 10: Enable and validate the Si64 adaptive block-CG route

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/run_dg_hybrid_si64_scf.py`
- Modify: `docs/verification/si64-hybrid-self-consistent-ground-state.md`

**Step 1: Add a failing comparison mode**

Require the block-CG run to consume the same fixed basis and initial DC+LCFO
density as the reference.  Parse adaptive tolerance, iteration count, stop
reason, Pulay reset, and memory receipts.

**Step 2: Connect block-CG behind an explicit solver option**

Default remains complex ScaLAPACK until the comparison passes.  Preserve the exact outer
SCF and checkpoint gates.

**Step 3: Run Si64 and compare physical results**

Compare converged density, occupied projector, total energy, electron count, and
symmetry receipts against Task 9.  Diagnose any oscillation using both inner and
outer histories; do not increase the CG cap to hide it.

**Step 4: Run the full focused suite**

```bash
python3 tests/dg/run_dg_hybrid_ground_state_types_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_density_mpi.py
python3 tests/dg/run_dg_hybrid_scf_mpi.py
python3 tests/dg/run_dg_hybrid_block_cg_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
git diff --check
```

Expected: PASS on MPI 1/2/4/8 where supported.

**Step 5: Independent review and commit**

Review numerical stability, collective safety, memory receipts, and the absence
of a hidden fixed high CG iteration count before committing.

```bash
git add src/gs/main_dft.f90 src/io/inputoutput.f90 \
  tests/dg/run_dg_hybrid_si64_scf.py \
  docs/verification/si64-hybrid-self-consistent-ground-state.md
git commit -m "feat: enable adaptive block CG for Si64 hybrid SCF"
```
