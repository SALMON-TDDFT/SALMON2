# Direct DG Term in Existing DC CG Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the projected local-basis DG ground-state route with the
existing DC orbital CG plus a six-face DG/SIPG Hamiltonian term.

**Architecture:** The conventional DC SCF loop remains the owner of orbitals,
orthogonalization, density, mixing, Hartree, XC, and local potential. The DG
operator is a linear add-on applied to both the current orbital and CG search
direction at every Hamiltonian evaluation. SCF-level continuation controls
only its scalar weight.

**Tech Stack:** Fortran 2008, MPI, OpenMP, CMake, Python contract fixtures.

---

### Task 1: Wire a direct DG operator into the existing DC orbital CG

**Files:**
- Modify: `src/gs/scf_iteration.f90`
- Modify: `src/gs/conjugate_gradient.f90`
- Modify: `src/gs/scf_iteration_dft.f90`
- Test: `tests/dg/check_dg_dc_local_periodic_route.py`
- Test: `tests/dg/test_dg_dc_direct_cg_mpi.f90`
- Create: `tests/dg/run_dg_dc_direct_cg_mpi.py`

**Step 1: Write RED contract and MPI tests**

Require:

- `solve_orbitals` passes `dc` and the DG scale only for
  `yn_dg_dc_local_periodic='y'`;
- `gscg_rwf` applies the DG action to `xk`, every `pk`, and every refreshed
  `xk`;
- all `system%no == nstate_frag` states are retained;
- the six-face operator is linear and zero at scale zero;
- no LCFO, EigenExa, projected metric, or global coefficient solve is called.

Run:

```bash
python3 tests/dg/check_dg_dc_local_periodic_route.py
python3 tests/dg/run_dg_dc_direct_cg_mpi.py
```

Expected: RED before wiring.

**Step 2: Reuse the existing face-action implementation**

Generalize the existing optional flux Hamiltonian hook in `gscg_rwf` into a
route-neutral additional interface action. Keep the `hpsi` volume/nonlocal
call unchanged, then add `lambda * H_DG psi`. Repeat for the search direction
and refreshed orbital. Do not allocate a global orbital or coefficient
matrix.

**Step 3: Keep standard DC density and potential ownership**

Leave these calls in their existing order:

```fortran
call calc_density(...)
call calc_rho_total_dcdft(...)
call update_density_and_potential(...)
call calc_vlocal_fragment_dcdft(...)
```

Do not call the local-basis density reconstruction or its projected
Hamiltonian routines from the production driver.

**Step 4: Run focused verification**

Run the two Task 1 tests, existing DC handoff MPI fixture, and
`git diff --check`. Expected: PASS.

**Step 5: Specification and code-quality review**

Resolve every Critical/Important finding before continuing.

### Task 2: Move continuation and checkpointing onto the direct DC state

**Files:**
- Modify: `src/gs/scf_iteration_dft.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/common/dg_ground_state_checkpoint.f90`
- Test: `tests/dg/test_dg_ground_state_checkpoint_mpi.f90`

**Step 1: Write RED tests**

Require an accepted scale-one unmixed fixed point, measured orbital/density,
orthogonality, Hermiticity, charge, and face diagnostics, and a checkpoint
containing fragment orbitals rather than projected coefficient rows.

**Step 2: Implement SCF-level continuation**

Advance the DG scale only between completed DC-style SCF stages. On a failed
stage restore orbitals, density, potentials, and scale. Recompute the DG
action during every inner CG Hamiltonian application.

**Step 3: Publish direct-state checkpoint**

Store fragment core-plus-buffer orbitals, occupations, density, Hartree, XC,
local potential, fingerprints, topology, continuation history, and measured
diagnostics. Keep standard checkpoint publication disabled.

**Step 4: Focused verification and two reviews**

Run checkpoint round-trip/corruption/rank-disagreement tests, route tests,
overlay build, specification review, and code-quality review. Resolve every
Critical/Important finding.

### Task 3: Re-run the Si64 gate

**Files:**
- Modify: `docs/plans/2026-07-24-dg-dc-local-periodic-wannier.md`

**Step 1: Establish the normal-DC baseline**

Use the verified Si64 input and `atom.dat`, `nstate_frag=400`,
`mixrate=0.2`, eight MPI ranks, and one OpenMP thread. Confirm the DC trace
matches the known 87-SCF `threshold=1e-9` run.

**Step 2: Run the handoff matrix**

Use fresh directories for `1e-2`, `1e-3`, and `1e-4`, at least two buffer
widths and two equivalent fragment decompositions. Include a fully converged
DC (`1e-9`) reference handoff.

**Step 3: Record all gates**

Record runtime, memory, state count, continuation/rollback history, orbital
and density residuals, charge, energy, orthogonality, Hermiticity, six-face
balance, and checkpoint publication.

**Step 4: Stop or commit**

If no Si64 configuration passes, document the failure and stop before
Wannier Tasks 6–8. Otherwise complete specification and code-quality reviews,
resolve Critical/Important findings, run a fresh parent-prerequisite overlay
build, and commit Task 5.
