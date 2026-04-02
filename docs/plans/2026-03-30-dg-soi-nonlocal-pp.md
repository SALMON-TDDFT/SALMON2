# DG-SOI Nonlocal PP Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add the full SALMON spin-orbit nonlocal pseudopotential contribution to the DG-SOI Hamiltonian and propagation path.

**Architecture:** Reuse the existing full-SOI pseudopotential data (`zekr_uV_so`, `rinv_uvu_so`) and build DG fragment-basis matrix elements directly in the SOI Hamiltonian builder. Keep scalar nonlocal PP and SO nonlocal PP as separate helpers and sum their contributions into the DG-SOI complex Hamiltonian path.

**Tech Stack:** Fortran, SALMON DG-RT modules, spin-orbit pseudopotential data prepared in `src/so`, MPI/OpenMP existing DG infrastructure

---

### Task 1: Audit DG-SOI Hamiltonian rebuild entry points

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90`
- Modify: `src/rt/dg/rt_dg_fragment_basis_update_soi.f90`
- Modify: `src/rt/dg/rt_dg_integrator_derivative.f90`

**Step 1: Identify all SOI Hamiltonian construction paths**

Confirm where DG-SOI Hamiltonian/nonlocal operator data are built or rebuilt:

- initial DG-SOI Hamiltonian construction
- adaptive basis update rebuild
- any cached nonlocal/operator refresh in time propagation

**Step 2: Add temporary comments marking required SO-NL insertion points**

Add short comments near each rebuild path noting that SO nonlocal PP must be applied there.

**Step 3: Build**

Run: `make -C build -j2 salmon`

Expected: PASS

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90 src/rt/dg/rt_dg_fragment_basis_update_soi.f90 src/rt/dg/rt_dg_integrator_derivative.f90
git commit -m "chore: mark DG-SOI SO nonlocal rebuild sites"
```

### Task 2: Add a DG-SOI SO nonlocal matrix helper

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90`

**Step 1: Write the helper skeleton**

Add a new routine, for example:

- `add_nonlocal_pp_matrix_so(...)`

It should take:

- `dg_frag`
- `mg`
- `ppg`
- `nspin`
- `hvol`

and immediately return if SO projector data are absent.

**Step 2: Build projector overlaps from `zekr_uV_so`**

Inside the helper, mirror the overlap-building structure used by scalar nonlocal PP, but use:

- `ppg%ia_tbl_so`
- `ppg%zekr_uV_so`
- `ppg%rinv_uvu_so`
- `phi_frag_c` for SOI basis overlaps

**Step 3: Accumulate complex matrix elements**

Accumulate the SO nonlocal contribution into the DG-SOI complex Hamiltonian storage used by the surrounding SOI path.

**Step 4: Build**

Run: `make -C build -j2 salmon`

Expected: PASS

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90
git commit -m "feat: add DG-SOI SO nonlocal PP matrix builder"
```

### Task 3: Wire SO nonlocal PP into initial DG-SOI Hamiltonian construction

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90`

**Step 1: Find the scalar nonlocal PP insertion point**

Locate where scalar nonlocal PP is added in the SOI Hamiltonian build.

**Step 2: Add SO nonlocal PP accumulation**

Call the new SO helper after the scalar nonlocal contribution so both are included in the final DG-SOI Hamiltonian.

**Step 3: Preserve Hermiticity/symmetry handling**

Ensure the existing post-processing still acts on the combined matrix.

**Step 4: Build**

Run: `make -C build -j2 salmon`

Expected: PASS

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90
git commit -m "feat: include SO nonlocal PP in DG-SOI Hamiltonian"
```

### Task 4: Wire SO nonlocal PP into SOI basis-update rebuilds

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_basis_update_soi.f90`

**Step 1: Find basis-update Hamiltonian rebuild path**

Locate the SOI basis-update routine that reconstructs or refreshes DG Hamiltonian operators.

**Step 2: Reuse the same SO nonlocal helper**

Ensure the basis-update path includes the same SO nonlocal PP contribution as initial construction.

**Step 3: Build**

Run: `make -C build -j2 salmon`

Expected: PASS

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_basis_update_soi.f90
git commit -m "fix: include SO nonlocal PP in DG-SOI basis rebuild"
```

### Task 5: Verify derivative/operator path consistency

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_derivative.f90`
- Modify: `src/rt/dg/rt_dg_fragment_ops.f90` (only if needed)

**Step 1: Check how the propagated Hamiltonian is sourced**

Confirm whether time evolution reads the rebuilt DG Hamiltonian directly or needs a separate cached nonlocal refresh.

**Step 2: Add the minimum consistency hook**

If no extra hook is needed, add a clarifying comment.

If a hook is needed, wire the SO nonlocal piece through the same cache/update path used by scalar nonlocal PP.

**Step 3: Build**

Run: `make -C build -j2 salmon`

Expected: PASS

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_integrator_derivative.f90 src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "fix: align DG-SOI derivative path with SO nonlocal PP"
```

### Task 6: Add lightweight diagnostics for first verification

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90`

**Step 1: Add guarded diagnostics**

Under a local debug flag or existing trace style, print:

- max absolute SO nonlocal matrix contribution
- whether SO projector arrays are allocated

for one root rank during initialization.

**Step 2: Build**

Run: `make -C build -j2 salmon`

Expected: PASS

**Step 3: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90
git commit -m "chore: add DG-SOI SO nonlocal verification trace"
```

### Task 7: Run a minimal SOI smoke verification

**Files:**
- Test: existing SOI DG-RT input case

**Step 1: Build final binary**

Run: `make -C build -j2 salmon`

Expected: PASS

**Step 2: Run a small SOI DG-RT case**

Run an existing small SOI DG-RT case with `yn_spinorbit='y'`.

Expected:

- no crash
- SO nonlocal diagnostic is nonzero for a pseudopotential with SOI
- Hamiltonian remains finite

**Step 3: Compare against previous behavior**

Check whether:

- SO nonlocal term is now reported
- energy/Hamiltonian shifts relative to the old DG-SOI path are plausible

**Step 4: Commit**

```bash
git add .
git commit -m "test: verify DG-SOI SO nonlocal PP integration"
```
