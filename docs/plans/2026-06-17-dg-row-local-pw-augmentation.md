# DG Row-Local PW Augmentation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Re-enable DG plane-wave helper channels without storing full mixed coefficient vectors or full fragment-PW operator rows on every rank.

**Architecture:** Add row-owner-local fragment-PW operator storage and a row-local mixed Hamiltonian apply. Keep the existing pure-fragment Taylor path intact, and only reconnect PW propagation after every full-gather path has been replaced by local row requests.

**Tech Stack:** Fortran 90/2003, MPI communication wrappers in `communication`, OpenMP where existing loops already use it, CMake build `build-mpi-eigenexa`.

---

### Task 1: Add Local PW Operator Fields

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment.f90`
- Modify: `src/rt/dg/rt_dg_fragment_soi.f90`
- Modify: `src/rt/dg/rt_dg_plane_wave.f90`

**Step 1: Add local-row fields to `s_dg_fragment_rt`**

Add fields near the existing PW arrays:

```fortran
integer, allocatable :: fp_local_row_ids(:)
integer, allocatable :: fp_local_pw_ids(:)
complex(8), allocatable :: S_mat_frag_pw_local(:,:,:)
complex(8), allocatable :: H_mat_frag_pw_local(:,:,:)
complex(8), allocatable :: P_mat_frag_pw_local(:,:,:,:)
```

Use dimensions:

```text
S/H: (n_local_fragment_rows, n_requested_pw_rows, nspin)
P:   (3, n_local_fragment_rows, n_requested_pw_rows, nspin)
```

**Step 2: Clear the fields when PW data is cleared**

In cleanup paths that already deallocate `S_mat_frag_pw` and `H_mat_frag_pw`, also deallocate the local arrays and row-id lists.

**Step 3: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_types.f90 src/rt/dg/rt_dg_fragment.f90 src/rt/dg/rt_dg_fragment_soi.f90 src/rt/dg/rt_dg_plane_wave.f90
git commit -m "feat: add row-local DG PW operator storage"
```

### Task 2: Populate Local Fragment-PW Blocks

**Files:**
- Modify: `src/rt/dg/rt_dg_plane_wave.f90`
- Modify: `src/rt/dg/rt_dg_integrator_stage_update.f90`

**Step 1: Add a helper to derive local fragment rows**

Create a helper in `rt_dg_plane_wave.f90`:

```fortran
subroutine build_local_fragment_pw_row_list(dg_frag, ispin)
```

It should fill `fp_local_row_ids` from `local_coef_global_ids` when available. If not available, use row-owner logic and include only rows owned by `dg_frag%id`.

**Step 2: Add a helper to derive requested PW rows**

For the first implementation, use all PW rows only if `n_plane_waves` is tiny and mark it as transitional. The allocation must still be local-row in the fragment dimension. Later tasks can distribute PW rows for `Y_P`.

**Step 3: Populate `S_mat_frag_pw_local`**

Reuse `compute_fragment_pw_overlap` logic but write only rows listed in `fp_local_row_ids`.

**Step 4: Populate `H_mat_frag_pw_local` and `P_mat_frag_pw_local`**

Reuse the existing mixed matrix formulas, but loop over local fragment rows only. Do not allocate or fill `n_mat_max x n_pw` arrays for RT propagation.

**Step 5: Wire stage update**

In `rt_dg_integrator_stage_update.f90`, replace the PW mixed update allocation for RT with local-row population. Keep dense/full mixed helpers only for DC diagnostics or disabled paths.

**Step 6: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 7: Commit**

```bash
git add src/rt/dg/rt_dg_plane_wave.f90 src/rt/dg/rt_dg_integrator_stage_update.f90
git commit -m "feat: build row-local fragment-PW blocks"
```

### Task 3: Implement Row-Local Mixed Hamiltonian Apply

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `src/rt/dg/rt_dg_integrator_derivative.f90`

**Step 1: Add a new operator apply**

Create:

```fortran
subroutine apply_mixed_hamiltonian_local_rows(dg_frag, ispin, coef_frag, coef_pw, y_frag, y_pw)
```

Inputs should be compact buffers, not full mixed vectors.

**Step 2: Fragment output contribution**

Compute:

```text
y_frag = H_FF(local, needed F) * coef_frag + H_FP(local, needed P) * coef_pw
```

Use existing compact `H_mat_blocks` apply for `H_FF`, and `H_mat_frag_pw_local` for `H_FP`.

**Step 3: PW output contribution**

For the first reconnect, support only PW rows owned by the rank. Accumulate `H_PF` using the Hermitian transpose of local `H_FP`, then reduce only owned PW rows if needed.

**Step 4: Replace derivative PW stop**

In `calculate_time_derivative`, replace:

```fortran
stop "DG derivative PW path requires row-local mixed apply; full-basis gather is disabled"
```

with the new local mixed derivative path.

**Step 5: Add all-row guard**

Keep a fatal guard if any PW derivative path requests `1:n_mat_max` fragment rows.

**Step 6: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 7: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_ops.f90 src/rt/dg/rt_dg_integrator_derivative.f90
git commit -m "feat: apply DG PW Hamiltonian on local rows"
```

### Task 4: Reconnect Taylor4-PC PW Propagation

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_taylor.f90`

**Step 1: Replace the PW stop**

Replace:

```fortran
stop "DG Taylor4-PC PW path requires row-local mixed apply; full-basis gather is disabled"
```

with a branch that allocates local fragment and local PW coefficient blocks.

**Step 2: Keep H fixed during Taylor**

Ensure stage update runs before the Taylor expansion and is not called inside the Taylor order loop.

**Step 3: Propagate only owned rows**

Update `dg_frag%coef` and `dg_frag%coef_pw` only for rows owned by this rank. Do not call `zero_nonowned_coefficients` as a replacement for correct local updates.

**Step 4: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_integrator_taylor.f90
git commit -m "feat: reconnect Taylor4-PC row-local PW propagation"
```

### Task 5: Replace PW Density Full Cache

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Modify: `src/rt/dg/rt_dg_fragment_ops.f90`

**Step 1: Remove the density PW stop**

Replace:

```fortran
stop "DG density PW path requires row-local PW reconstruction; full PW coefficient cache is disabled"
```

with a row-local density path.

**Step 2: Add batched PW coefficient fetch**

Fetch only PW rows needed for the current density batch. Do not repopulate `coef_pw_full_cache`.

**Step 3: Accumulate density by grid block**

For each grid block, compute PW wavefunction contribution from fetched PW rows and reduce density, not coefficients.

**Step 4: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90 src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "feat: reconstruct DG PW density without full coefficient cache"
```

### Task 6: Re-enable Diagnostics Safely

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_unitarity.f90`
- Modify: `src/rt/dg/rt_dg_observables.f90`

**Step 1: Replace unitarity PW stop**

Use row-local overlap apply. If a global norm is needed, reduce scalar overlaps, not vectors.

**Step 2: Audit observables**

Search:

```bash
rg -n "coef_pw_full_cache|gather_full_coef_view|1:n_mat_max|n_mat_max, n_plane_waves" src/rt/dg
```

Expected: no active RT path gathers full mixed state.

**Step 3: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_integrator_unitarity.f90 src/rt/dg/rt_dg_observables.f90
git commit -m "feat: make DG PW diagnostics row-local"
```

### Task 7: Smoke and Weak-Scaling Checks

**Files:**
- Modify or create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_pw_smoke`
- Optional create: `docs/plans/2026-06-17-dg-row-local-pw-validation.md`

**Step 1: Run static checks**

Run:

```bash
git diff --check
```

Expected: no output.

**Step 2: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 3: Run tiny PW smoke**

Run a short RT case with small `n_plane_waves_dg`, for example 4 or 8.

Expected:

- no full-gather fatal stops,
- no NaN in density or current,
- output includes finite `J_tot,z`.

**Step 4: Compare 2x2x2 and 8x8x8 behavior**

Use the same per-fragment basis policy and measure timing sections.

Expected:

- no obvious `n_frag_total` scaling in Taylor derivative,
- no replicated PW coefficient cache allocation,
- current behavior is at least physically inspectable.

**Step 5: Commit validation notes**

```bash
git add samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_pw_smoke docs/plans/2026-06-17-dg-row-local-pw-validation.md
git commit -m "test: validate row-local DG PW smoke path"
```
