# Neighbor-Coupled Fragment-Local RT Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement and validate a correctness-first fragment-local RT backend that propagates each owner fragment together with face-neighbor cross couplings, then writes back only the owner component.

**Architecture:** Add a namelist-selected diagnostic backend beside the existing global mixed-Z and isolated fragment-local routes. Reuse the existing face-neighbor enumeration and dense `expdiag` machinery where possible. Build an environment block containing owner plus face-neighbor coefficients, assemble `S/H/R`, propagate the environment block, and owner-writeback the result.

**Tech Stack:** Fortran 90 in `src/rt/dg`, SALMON namelist plumbing in `src/io`, existing DG fragment/BPW data structures, MPI coefficient synchronization, dense diagonalization through the current expdiag/eigen solver path, Python/NumPy analysis scripts for Pz and HHG comparisons.

---

### Task 1: Audit Existing Neighbor And Matrix Data

**Files:**
- Inspect: `src/rt/dg/rt_dg_wpw_window.f90`
- Inspect: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Inspect: `src/rt/dg/rt_dg_plane_wave.f90`
- Inspect: `src/rt/dg/rt_dg_mixed_fsum_diagnose.f90`
- Inspect: `src/rt/dg/rt_dg_fragment*.f90`

**Step 1: Locate the canonical face-neighbor helper**

Confirm that the new backend can reuse `wpw_face_neighbor_fragment` or the equivalent runtime helper, and document any mismatch with Hamiltonian-side neighbor enumeration.

**Step 2: Locate available owner-neighbor blocks**

Check which of the following are already available at RT initialization:

- overlap blocks `S_owner,neighbor`
- field-free Hamiltonian blocks `H_owner,neighbor`
- position blocks `R_owner,neighbor`
- coefficient storage for owner and neighbor slots

**Step 3: Add a read-only diagnostic if needed**

If the available data path is unclear, add a temporary or gated diagnostic routine that reports whether each cross block exists, its dimensions, and Hermiticity against the transpose-conjugate partner.

**Step 4: Build**

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: build succeeds.

### Task 2: Add A Namelist Backend Selector

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: the DG RT initialization file that selects `dg_mixed_z_local_prop_backend`

**Step 1: Add a static test**

Check that the namelist contains a value such as:

```text
dg_mixed_z_local_prop_backend = 'neighbor_env_expdiag'
```

Expected before implementation: FAIL.

**Step 2: Wire the selector**

Add the new backend string behind existing DG/mixed-Z controls. Do not add an environment-variable-only switch.

**Step 3: Route to a stub backend**

The stub should stop with a clear message if called before Task 4 is implemented:

```text
[DG-MIXED-Z-NEIGHBOR-ENV] backend selected but propagation kernel is not implemented
```

**Step 4: Build**

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: build succeeds.

### Task 3: Build Environment Layout And Coefficient Synchronization

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Modify or add helper under: `src/rt/dg`

**Step 1: Add layout diagnostics**

For each owner fragment, print in diagnostic mode:

- owner fragment ID
- face-neighbor IDs
- owner slot count
- neighbor slot count
- total environment dimension

**Step 2: Gather coefficients every step**

Create the environment coefficient block from the current owner and face-neighbor coefficients. The gather must occur every RT step before propagation so that neighbor states are current.

**Step 3: Owner-only writeback**

After environment propagation, copy only the owner component back to the global owned coefficient storage.

**Step 4: Field-off smoke**

Run a short field-off calculation and confirm:

- no fatal missing-neighbor errors
- coefficient norms remain finite
- owner-only writeback does not double-count neighbor slots

### Task 4: Assemble `S_env` And `H_env`

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Modify or add helper under: `src/rt/dg`

**Step 1: Assemble block matrices**

Build:

```text
S_env = [[S_oo, S_on],
         [S_no, S_nn]]

H_env = [[H_oo, H_on],
         [H_no, H_nn]]
```

Use the same basis ordering as the environment coefficient block.

**Step 2: Add Hermiticity and metric diagnostics**

Report:

- `max_abs(S - S^H)`
- `max_abs(H - H^H)`
- smallest metric eigenvalue or equivalent stability diagnostic

**Step 3: Generalized expdiag propagation**

Use the existing dense generalized diagonalization path if available. If the current expdiag route assumes an orthonormal basis, add a local orthonormalization step for the environment block.

**Step 4: Field-off comparison**

Compare against:

- isolated fragment-local route
- global mixed-Z reference

Expected: field-off drift and coefficient mismatch should improve relative to isolated local propagation or at least not worsen.

### Task 5: Assemble Field-Coupling `R_env`

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Modify: the mixed-Z position reconstruction/export path if cross `R` blocks are absent

**Step 1: Confirm cross position availability**

If `R_owner,neighbor` and `R_neighbor,owner` are missing, add the minimal reconstruction/export path needed for face-neighbor environment propagation.

**Step 2: Build the field-dependent Hamiltonian**

For length-gauge propagation use:

```text
H_eff(t) = H_env - E_z(t) Rz_env
```

Keep the impulse-specific momentum-kick route separate from this laser route.

**Step 3: Add diagnostics**

Report `R_env` Hermiticity and compare the owner-only `R_oo` result against the full environment `R_env` result for a short run.

### Task 6: Validate C64 Laser Odd Response

**Files:**
- Use C64 fragment-local sample inputs under `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac`
- Add temporary copied inputs only if needed

**Step 1: Run short `+E` and `-E` laser tests**

Use the same C64 settings that previously showed large even-order peaks. Keep `OMP_NUM_THREADS=2`.

**Step 2: Plot time-domain `Pz`**

Check that:

- `Pz(+E) + Pz(-E)` is small
- the offset after the pulse is reduced
- the waveform no longer oscillates once per half-cycle of `Az`

**Step 3: Compute HHG**

Use the same FFT normalization and windowing as the current comparison scripts. Report:

- `H1`
- `H2/H1`
- `H3/H1`
- `H4/H1`
- `H5/H1`

Expected: `H2` and `H4` should become dips for centrosymmetric C64.

### Task 7: Compare Against Global Reference And Decide Next Route

**Files:**
- Use existing global mixed-Z reference inputs/outputs where possible
- Add a plot under `/private/tmp` or a tracked docs artifact only if it is intentionally part of the evidence

**Step 1: Overlay local-neighbor and global `Pz`**

Compare amplitude, phase, offset, and odd-response symmetry.

**Step 2: Overlay HHG spectra**

The neighbor-coupled route does not need to be fast yet, but it should reproduce the physically important low-order structure better than isolated fragment-local propagation.

**Step 3: Commit the implementation**

Commit only the source/input changes required for the backend and the small evidence scripts or docs. Leave generated run outputs untracked unless the user explicitly asks to archive them.

**Step 4: Decide optimization**

If the dense environment route is physically correct, plan a later sparse/weak-scaling implementation. If it still fails odd symmetry, inspect the exported cross-position blocks and the owner-writeback algebra before changing Wannier generation again.
