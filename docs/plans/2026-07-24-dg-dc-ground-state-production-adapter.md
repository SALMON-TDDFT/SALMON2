# DG-DC Ground-State Production Adapter Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Complete Task 4 by invoking adaptive DG continuation from `main_dft` with a GS-native operator and the existing DC density/Hartree/LDA-XC/vlocal route.

**Architecture:** A DC-owned adapter maps common nodal orbitals to existing fragment and total-system arrays. It applies the always-full volume/nonlocal Hamiltonian plus a lambda-scaled physical SIPG face action, while density mixing remains transactional and the established DC potential routines remain authoritative.

**Tech Stack:** Fortran 2008, MPI, SALMON DC/Poisson/XC/pseudopotential infrastructure, common nodal SIPG API, Python source contracts, CMake.

---

### Task A: Freeze the production adapter contract

**Files:**
- Create: `tests/dg/test_dg_dc_ground_state_adapter_mpi.f90`
- Create: `tests/dg/run_dg_dc_ground_state_adapter_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`
- Modify: `src/common/dg_nodal_interfaces.f90`

**Step 1: Write RED tests**

Require an adapter callback that receives independent volume/nonlocal and
physical-SIPG actions, candidate validity/occupations, and a DC potential
update callback. Verify that only SIPG is lambda-scaled, padded candidates are
excluded, and density/potential reset is transactional.

**Step 2: Run RED**

```bash
python3 tests/dg/run_dg_dc_ground_state_adapter_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
```

Expected: FAIL because the production adapter and `main_dft` invocation do not
exist.

**Step 3: Add minimal common callback contracts**

Expose separate volume/nonlocal and SIPG actions. Keep SALMON/DC structures out
of `src/common`.

**Step 4: Verify the fixture compiles to the intended missing symbol**

Re-run the RED commands and record the failure.

### Task B: Implement the GS-native production callback

**Files:**
- Create: `src/gs/dc/dg_dc_ground_state_adapter.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `tests/dg/test_dg_dc_ground_state_adapter_mpi.f90`

**Step 1: Implement Hamiltonian composition**

Compute:

```fortran
hpsi = hpsi_volume_nonlocal + lambda*hpsi_sipg
```

and derive orbital, Hermiticity, orthogonality, physical-face, and charge
diagnostics from the same epoch.

**Step 2: Implement candidate ownership**

Use explicit valid-candidate masks and occupations. Exclude MPI padding from
orthogonalization, projector tracking, and density.

**Step 3: Implement density and potential callback order**

Reconstruct core-owned fragment density, linearly mix against the accepted
density, call the bound DC density/potential update, and retain both old and
new potential for rollback.

**Step 4: Run focused adapter tests**

```bash
python3 tests/dg/run_dg_dc_ground_state_adapter_mpi.py
python3 tests/dg/run_dg_dc_ground_state_mpi.py
python3 tests/dg/run_dg_nodal_sipg_mpi.py
```

Expected: PASS.

### Task C: Bind existing DC density, Hartree, XC, and vlocal

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/scf_iteration.f90`
- Modify: `src/gs/dc/dcdft.f90`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Expose the existing DC potential update safely**

Refactor the reusable body of `update_density_and_potential` behind an
explicit interface without changing the normal SCF caller or numerical order.

**Step 2: Add the `main_dft` host binding**

Immediately after `materialize_dg_dc_candidates_for_main`:

1. map nodal core density into fragment `rho_s`;
2. call `calc_rho_total_dcdft`;
3. apply the configured simple linear mix;
4. call the existing Hartree/LDA-XC/local-potential update;
5. call `calc_vlocal_fragment_dcdft`;
6. invoke `run_dg_dc_ground_state` collectively on `dc%icomm_tot`.

Fail collectively unless the result is accepted, at lambda one, and unmixed
verified. Do not invoke LCFO/WPW/checkpoint/RT.

**Step 3: Run focused route tests**

```bash
python3 tests/dg/check_dg_dc_local_periodic_route.py
python3 tests/dg/run_dg_dc_ground_state_adapter_mpi.py
python3 tests/dg/run_dg_dc_ground_state_mpi.py
```

Expected: PASS.

### Task D: Complete Task 4 review and verification

**Files:**
- Review all Task 4 and adapter changes.

**Step 1: Run required regressions**

Run continuation, SIPG, canonical-face, complete-H, nonlocal, density, and
normal DC route fixtures, including bounds checking where provided.

**Step 2: Request specification review**

Review production invocation, lambda mathematics, DC potential reuse,
candidate ownership, rollback, MPI order, and non-promotion constraints.
Resolve every Critical/Important finding.

**Step 3: Request code-quality review**

Review dependency direction, error handling, numerical stability, state
ownership, cleanup, and maintainability. Resolve every Critical/Important
finding.

**Step 4: Build the parent-prerequisite overlay**

Copy the latest Task 4 files into a fresh `/tmp` parent-prerequisite overlay
without importing parent uncommitted changes into the branch. Configure and
build Release with MPI, ScaLAPACK, and EigenExa.

Expected: `[100%] Built target salmon`.

**Step 5: Verify and commit**

```bash
git diff --check
git commit -m "feat(dg): add self-consistent DC-to-DG continuation"
```

Commit only after all required verification and reviews pass.
