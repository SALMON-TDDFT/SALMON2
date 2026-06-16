# DG PW-Augmented Basis Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the abandoned real-space RT direction with a production-ready coefficient-space DG path that augments DC-LCFO-Flux bases using plane-wave helper functions.

**Architecture:** Keep RT propagation in DG coefficient space. First revert the real-space RT commits, then audit and enable the existing PW augmentation hooks by making overlap, Hamiltonian, derivative/Taylor, density, and current paths use the same generalized mixed basis consistently.

**Tech Stack:** Fortran 2008, SALMON DG-RT modules, MPI, OpenMP, EigenExa-enabled build, existing `rt_dg_plane_wave.f90` helpers.

---

### Task 1: Revert The Abandoned Real-Space RT Path

**Files:**
- Revert commits:
  - `7e5a5ab feat: add complex DC LCFO flux operator`
  - `1bd8dde refactor: share DC LCFO flux operator`
  - `269b9ca feat: add DG flux real-space RT mode flag`
  - `3309019 docs: plan DG flux real-space RT implementation`
  - `951f7a0 docs: design DG flux real-space RT path`

**Step 1: Confirm clean tracked state**

Run:

```bash
git status --short
```

Expected: only untracked build/sample directories, no tracked modifications.

**Step 2: Revert the commits newest-first**

Run:

```bash
git revert --no-edit 7e5a5ab
git revert --no-edit 1bd8dde
git revert --no-edit 269b9ca
git revert --no-edit 3309019
git revert --no-edit 951f7a0
```

Expected: each revert completes without conflict.

**Step 3: Verify removed symbols**

Run:

```bash
rg -n "yn_dg_flux_realspace_rt|dc_lcfo_flux_operator|apply_dc_lcfo_flux_hpsi_zwf|DG Flux Real-Space RT" src docs/plans
```

Expected: no matches, except historical text in git log is not searched.

**Step 4: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

---

### Task 2: Audit Existing PW Augmentation Wiring

**Files:**
- Inspect: `src/rt/dg/rt_dg_plane_wave.f90`
- Inspect: `src/rt/dg/rt_dg_fragment.f90`
- Inspect: `src/rt/dg/rt_dg_fragment_ops.f90`
- Inspect: `src/rt/dg/rt_dg_integrator_derivative.f90`
- Inspect: `src/rt/dg/rt_dg_integrator_taylor.f90`
- Inspect: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Inspect: `src/rt/dg/rt_dg_observables.f90`

**Step 1: Locate all PW stop gates**

Run:

```bash
rg -n "use_plane_wave_basis|coef_pw|PW|plane wave|supports the pure fragment" src/rt/dg
```

Expected: list all branches that already handle PW and all branches that currently stop.

**Step 2: Create an audit note**

Create `docs/plans/2026-06-17-dg-pw-augmented-basis-audit.md` with sections:

```markdown
# DG PW-Augmented Basis Audit

## Existing Complete Paths

## Stop Gates

## Missing Matrix Blocks

## Validation Inputs
```

**Step 3: Commit audit note**

```bash
git add docs/plans/2026-06-17-dg-pw-augmented-basis-audit.md
git commit -m "docs: audit DG PW augmentation wiring"
```

---

### Task 3: Add A Minimal PW Smoke Input

**Files:**
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_pw_smoke`

**Step 1: Copy existing short RT input**

Start from:

```bash
cp samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_smoke \
   samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_pw_smoke
```

**Step 2: Enable a small PW basis**

Set in `&propagation`:

```text
yn_plane_wave_basis = 'y'
n_plane_waves_dg = 8
k_cutoff_plane_wave = 0.1
yn_dg_subspace_diag = 'y'
dg_subspace_extra_states = 8
dg_subspace_pw_vectors = 4
```

Keep `nt` very short, for example `nt=2` or `nt=5`.

**Step 3: Verify RED**

Run:

```bash
mpirun -np 8 build-mpi-eigenexa/salmon < inputfile_rt_pw_smoke > /tmp/dg_pw_smoke_red.log 2>&1
```

Expected at this stage: fail at one of the existing PW stop gates, likely in `rt_dg_integrator_derivative.f90` or `rt_dg_integrator_taylor.f90`.

**Step 4: Commit smoke input**

```bash
git add samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_pw_smoke
git commit -m "test: add DG PW smoke input"
```

---

### Task 4: Make Taylor/Derivative Use Mixed Basis Apply

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_derivative.f90`
- Modify: `src/rt/dg/rt_dg_integrator_taylor.f90`
- Inspect/possibly modify: `src/rt/dg/rt_dg_plane_wave.f90`
- Inspect/possibly modify: `src/rt/dg/rt_dg_fragment_ops.f90`

**Step 1: Identify the mixed-basis apply helper**

Search:

```bash
rg -n "mixed|H_frag_pw|H_mat_pw|apply_.*pw|coef_pw|assemble_mixed" src/rt/dg
```

Expected: locate existing helper routines for mixed Hamiltonian application.

**Step 2: Replace PW stop in derivative path**

In `calculate_time_derivative`, replace:

```fortran
if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) then
  stop "DG derivative supports the pure fragment block-sparse route only"
end if
```

with a call to the mixed-basis derivative path. The minimum acceptable implementation should:

- apply `H_DC,DC`, `H_DC,PW`, `H_PW,DC`, and `H_PW,PW`
- include vector-potential momentum/current terms consistently
- use overlap/generalized basis convention already used by PW helpers

**Step 3: Replace PW stop in Taylor4-PC path**

In `time_evolution_taylor4pc`, replace:

```fortran
if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) then
  stop "DG Taylor4-PC supports the pure fragment block-sparse route only"
end if
```

with a mixed-basis Taylor path. Keep the Hamiltonian fixed for the Taylor operation.

**Step 4: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 5: Run PW smoke**

Run:

```bash
mpirun -np 8 build-mpi-eigenexa/salmon < samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_pw_smoke > /tmp/dg_pw_smoke_green.log 2>&1
```

Expected: no PW stop gate; run reaches at least `itt=2`.

**Step 6: Commit**

```bash
git add src/rt/dg/rt_dg_integrator_derivative.f90 src/rt/dg/rt_dg_integrator_taylor.f90 src/rt/dg/rt_dg_plane_wave.f90 src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "feat: propagate DG PW-augmented coefficients"
```

---

### Task 5: Verify Density And Current Consistency For PW Coefficients

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Modify: `src/rt/dg/rt_dg_observables.f90`
- Modify if needed: `src/rt/dg/rt_dg_fragment_ops.f90`

**Step 1: Run PW smoke with traces**

Run:

```bash
SALMON_DG_PW_TRACE=1 SALMON_DG_TRACE_CURRENT=1 \
mpirun -np 8 build-mpi-eigenexa/salmon < samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_pw_smoke > /tmp/dg_pw_density_current_red.log 2>&1
```

Expected: identify whether density and current include `coef_pw` and mixed terms.

**Step 2: Fix density reconstruction if incomplete**

Ensure density uses:

```text
rho = |psi_DC + psi_PW|^2
```

including:

- DC-DC
- DC-PW cross terms
- PW-PW terms

Do not compute redundant cross terms twice unless intentionally using symmetry.

**Step 3: Fix current observables if incomplete**

Ensure current/momentum uses the same mixed-basis coefficients:

- DC-DC momentum blocks
- DC-PW blocks
- PW-PW blocks
- nonlocal current contribution if available

**Step 4: Build and smoke**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
mpirun -np 8 build-mpi-eigenexa/salmon < samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_pw_smoke > /tmp/dg_pw_density_current_green.log 2>&1
```

Expected: run completes the short RT and writes current data.

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90 src/rt/dg/rt_dg_observables.f90 src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "fix: include PW coefficients in DG density and current"
```

---

### Task 6: Add PW Basis Diagnostics

**Files:**
- Modify: `src/rt/dg/rt_dg_plane_wave.f90`
- Modify: `src/rt/dg/rt_dg_fragment.f90`

**Step 1: Add root-rank diagnostics**

Print after PW basis setup:

```text
[DG-PW] requested=<n_plane_waves_dg> selected=<n_selected> kept=<n_plane_waves>
[DG-PW] k_cutoff=<...> highest_energy=<...>
[DG-PW] overlap_norm_dc_pw=<...> overlap_norm_pw_pw=<...>
```

**Step 2: Add condition diagnostics if overlap matrix is built**

Print min/max diagonal and a simple norm ratio. If a full eigenvalue condition estimate is already available, print min/max eigenvalue.

**Step 3: Build**

Run:

```bash
cmake --build build-mpi-eigenexa -j 8
```

Expected: build succeeds.

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_plane_wave.f90 src/rt/dg/rt_dg_fragment.f90
git commit -m "diag: report DG PW basis diagnostics"
```

---

### Task 7: Run Small PW Basis Sweep

**Files:**
- Test only unless failures reveal a code defect.

**Step 1: Prepare temporary inputs**

Create `/tmp/dg_pw_sweep` with copies of the smoke input using:

```text
n_plane_waves_dg = 0, 8, 16, 32
k_cutoff_plane_wave = 0.1
```

For `0`, use `yn_plane_wave_basis='n'`.

**Step 2: Run all cases**

Run each with:

```bash
mpirun -np 8 build-mpi-eigenexa/salmon < inputfile > run.log 2>&1
```

**Step 3: Extract current summary**

For each run, collect:

```text
J_tot,z at itt=1/2
J_para,z if trace is enabled
J_dia,z if trace is enabled
max |J_tot,x/y|
PW selected/kept count
```

**Step 4: Report**

Expected success trend:

- increasing PW count should increase paramagnetic cancellation
- `J_tot,z` should move closer to zero-centered response
- transverse current should not grow sharply

Do not claim final physics correctness from the smoke sweep alone.

---

### Task 8: Longer Diamond Impulse Check

**Files:**
- Test only unless failures reveal a code defect.

**Step 1: Choose best short-sweep PW size**

Pick the smallest PW count that improves cancellation without obvious instability.

**Step 2: Run longer impulse**

Use the established long RT input length from the previous diamond tests.

**Step 3: Assess qualitative physics**

Check:

- initial `J_tot,z`
- whether `J_tot,z` drops toward zero
- max/min symmetry
- whether oscillation center is near zero rather than `J_dia`
- transverse leakage

**Step 4: Report**

Summarize whether PW augmentation is addressing the f-sum/paramagnetic response problem.
