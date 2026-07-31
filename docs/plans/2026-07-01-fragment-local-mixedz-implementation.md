# Fragment-Local Mixed-Z Propagation Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement and validate a `W + Pself + Pneighbor` fragment-local mixed-Z propagation backend for DG real-time TDDFT.

**Architecture:** Keep `global_mixed_split_backend` as the reference implementation. Extend `fragment_local_mixed_split_backend` so it builds local mixed blocks, applies the field kick and phase locally, writes back owner coefficients, and logs coefficient/observable differences against the global backend until it is trusted.

**Tech Stack:** Fortran 90/95 in SALMON RT-DG modules, MPI reductions through existing `comm_summation`, EigenExa/LAPACK wrapper `eigen_zheev`, existing Si64 temporary benchmark inputs under `/private/tmp/stage2c_si64_impulse_pseudochannels_fixed`.

---

### Task 1: Baseline Diagnostics For Existing Local Backend

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Test input: `/private/tmp/stage2c_si64_impulse_pseudochannels_fixed/inputfile_si64_dg_impulse_z_dt01_nt4_fieldkickdiag`

**Step 1: Run the current local backend in diagnostic mode**

Run:

```bash
cd /private/tmp/stage2c_si64_impulse_pseudochannels_fixed
OMP_NUM_THREADS=2 \
SALMON_DG_EXPDIAG_GLOBAL_FLUX=1 \
SALMON_DG_EXPDIAG_GLOBAL_FIELD=1 \
SALMON_DG_MIXED_FSUM=1 \
mpiexec -n 8 /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/build-mpi-eigenexa-wannier-lib/salmon \
  < inputfile_si64_dg_impulse_z_dt01_nt4_fieldkickdiag \
  2>&1 | tee run_local_backend_baseline.log
```

Expected: global backend still runs. Existing local backend logs either unavailable or diagnostic-only status.

**Step 2: Add one concise local-backend summary line**

In `apply_fragment_local_mixed_split_exp_stub`, ensure every call prints one line containing:

- `step`
- `field_block_kind`
- `field_abs_sum`
- `available`
- `replacement_applied`
- `bad`
- `replacement_block_reason`
- `W_owner_storage_slots`
- `P_self_storage_slots`
- `P_neighbor_storage_slots`

**Step 3: Build**

Run:

```bash
cd /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: build succeeds.

### Task 2: Make Local Block Construction Explicit

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Inspect: `src/rt/dg/rt_dg_fragment_types.f90`

**Step 1: Factor local slot enumeration**

Create or refactor an internal helper near the existing local backend code:

```fortran
subroutine enumerate_fragment_local_mixed_slots(w_slot_count, pself_slot_count, pneighbor_slot_count, bad)
```

It should use the existing face-neighbor logic and count:

- owned W slots
- P slots on the owner fragment
- P slots on face-neighbor fragments

Expected: counts match the existing diagnostic values.

**Step 2: Add global ID arrays for local blocks if missing**

Ensure W and P storage records the global mixed-Z ID used by the reference backend:

- W global ID: `gid`
- P global ID: `nw + gp`

Expected: local blocks can be compared directly to `cmix_ref(global_id,state,spin)`.

**Step 3: Build**

Run the same build command as Task 1.

Expected: build succeeds.

### Task 3: Implement Field-Off Local Phase For W + Pself + Pneighbor

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`

**Step 1: Write the diagnostic comparison first**

In the local backend, compare the local field-off result against `apply_global_mixed_split_exp` for `E_use = 0`.

Expected metrics:

- `coef_local_vs_global_Snorm`
- `global_ref_Snorm`
- relative difference
- separate W, Pself, Pneighbor differences

**Step 2: Implement minimal field-off phase**

For each local slot with global mixed ID `g`, apply:

```fortran
phase_c = cos(dg_frag%mixed_wannier_bpw_eval(g, ispin) * dt)
phase_s = sin(dg_frag%mixed_wannier_bpw_eval(g, ispin) * dt)
c_local = cmplx(phase_c, -phase_s, kind=8) * c_local
```

This should match the global diagonal phase path.

**Step 3: Verify with no-field short input**

Create a temporary no-field `nt=4` input from the Si64 impulse input by setting the field to zero or using an existing no-field input if present.

Run with:

```bash
dg_mixed_z_local_prop_backend = 'fragment_local_mixed_split_backend'
dg_mixed_z_frag_local_field_block = 'all'
```

Expected: local-vs-global field-off coefficient difference is near numerical roundoff.

### Task 4: Implement Field-On Local Position Kick

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`

**Step 1: Build local `Z_local` matrix**

For each fragment-local block, construct:

```fortran
field_h(i,j) =
  - E_use(1) * mixed_wannier_bpw_z(1, gid_i, gid_j, ispin)
  - E_use(2) * mixed_wannier_bpw_z(2, gid_i, gid_j, ispin)
  - E_use(3) * mixed_wannier_bpw_z(3, gid_i, gid_j, ispin)
```

Use the local global-ID arrays for `gid_i` and `gid_j`.

**Step 2: Exponentiate the local field block**

Use the same `eigen_zheev` pattern as the global backend:

```fortran
call eigen_zheev(field_vec, field_eval, field_h)
tmp = matmul(conjg(transpose(field_h)), c_block)
tmp(i,:) = exp(-i * field_eval(i) * dt) * tmp(i,:)
c_block = matmul(field_h, tmp)
```

**Step 3: Apply field-free phase after the kick**

Apply the same diagonal phase as Task 3 after the field kick, matching the current split order.

**Step 4: Verify with `nt=4` impulse**

Run:

```bash
cd /private/tmp/stage2c_si64_impulse_pseudochannels_fixed
OMP_NUM_THREADS=2 \
SALMON_DG_EXPDIAG_GLOBAL_FLUX=1 \
SALMON_DG_EXPDIAG_GLOBAL_FIELD=1 \
SALMON_DG_MIXED_FSUM=1 \
mpiexec -n 8 /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/build-mpi-eigenexa-wannier-lib/salmon \
  < inputfile_si64_dg_impulse_z_dt01_nt4_fieldkickdiag
```

Expected: field-on local-vs-global differences are finite, logged, and small enough to guide the next fix. Do not switch production yet unless the difference is already acceptable.

### Task 5: Write Back Owner Coefficients Without Double Counting

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`

**Step 1: Write back W owner coefficients**

Use existing `writeback_fragment_local_storage_fieldoff` as the starting point. Ensure W slots are summed over owner fragments only and reconstructed through `global_wannier_coef`.

**Step 2: Write back Pself coefficients**

For each P global ID, write only owner `Pself` contributions to `mixed_wannier_bpw_pcoef`.

**Step 3: Do not write back Pneighbor as owner data**

Use Pneighbor for local field coupling, but write back only the owning fragment's P slot to avoid double counting.

**Step 4: Verify coefficient norm conservation**

Run `nt=4` impulse and inspect:

- norm before/after
- unoccupied weight
- local-vs-global coefficient difference

Expected: norm conservation comparable to global backend.

### Task 6: Compare Observables Against Global Backend

**Files:**
- Modify: `src/rt/dg/rt_dg_observables.f90` only if additional logs are needed
- Use existing: `src/rt/dg/rt_dg_fragment_ops.f90`

**Step 1: Run global backend reference**

Use corrected mixed-P input:

```bash
inputfile_si64_dg_impulse_z_dt01_nt2000_pw128_reseed_mixedP
```

**Step 2: Run local backend candidate**

Create a local-backend copy with:

```fortran
dg_mixed_z_local_prop_backend = 'fragment_local_mixed_split_backend'
dg_mixed_z_frag_local_field_block = 'all'
```

**Step 3: Compare Pz time series**

Use Python:

```python
import numpy as np
g = np.loadtxt("global_dg_polarization.data", comments="#")
l = np.loadtxt("local_dg_polarization.data", comments="#")
diff = l[:,3] - g[:,3]
print(np.sqrt(np.mean(diff**2)), np.max(np.abs(diff)))
```

Expected: RMS and max difference are small relative to the Pz signal. If not, inspect W/Pself/Pneighbor diagnostic differences.

### Task 7: Compare Si Impulse Spectra

**Files:**
- Use: `*_dg_response_from_p.data`

**Step 1: Run local backend `nt=2000` impulse**

Expected: `end SALMON`.

**Step 2: Compare `eps_z`**

Compare columns:

- Re eps z: column 10
- Im eps z: column 13

Windows:

- 0-1 eV
- 1-3 eV
- 3-5 eV
- 5-10 eV

Expected: local backend preserves the corrected DG mixed-P absorption structure.

### Task 8: Compare Si Laser HHG

**Files:**
- Use: `*_dg_polarization.data`

**Step 1: Run local backend laser**

Use the Si64 laser input:

```bash
inputfile_si64_dg_laser_z_dt01_nt4500_pw128_mixedP
```

with local backend selected.

**Step 2: Compute SALMON-like HHG from Pz**

Use:

```text
|J_z(omega)|^2 = omega^2 |P_z(omega)|^2
```

Expected: H1 remains dominant. H3/H1 and H7/H1 are comparable to the global backend order of magnitude. Large new even-order peaks indicate local block/writeback asymmetry.

### Task 9: Make Backend Selection A Namelist-Controlled Route

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`

**Step 1: Ensure existing namelist key is authoritative**

Use:

```fortran
dg_mixed_z_local_prop_backend = 'global_mixed_split_backend'
dg_mixed_z_local_prop_backend = 'fragment_local_mixed_split_backend'
```

Expected: no environment variable is required for normal backend selection.

**Step 2: Keep environment variables diagnostic-only**

Do not require `SALMON_DG_*` variables except for optional trace output.

### Task 10: Final Verification

**Files:**
- No code changes

**Step 1: Build**

Run:

```bash
cd /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: build succeeds.

**Step 2: Run smoke tests**

Run:

```bash
cd /private/tmp/stage2c_si64_impulse_pseudochannels_fixed
OMP_NUM_THREADS=2 mpiexec -n 8 /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/build-mpi-eigenexa-wannier-lib/salmon \
  < inputfile_si64_dg_impulse_z_dt01_nt4_fieldkickdiag
```

Expected: `end SALMON`.

**Step 3: Run local-backend `nt=4` impulse**

Expected: `end SALMON`, no nonfinite coefficients, logged local-vs-global differences.

**Step 4: Run process check**

Run:

```bash
ps -axo pid,command | rg -i "salmon|mpiexec|hydra" | rg -v "rg -i|Codex"
```

Expected: no remaining SALMON/MPI processes.
