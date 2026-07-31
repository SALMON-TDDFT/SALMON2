# Owner-Local Mixed-Z Redefinition Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the fragment-local mixed-Z production route with an owner-unique `W_owner + P_owned` local propagation space.

**Architecture:** The production backend should propagate only coefficients that have a unique owner and a unique writeback destination. Neighbor P channels can remain in diagnostics and future matrix-embedding experiments, but they must not be propagated and discarded in the production backend.

**Tech Stack:** Fortran 2008, SALMON DG fragment RT, MPI, existing CMake build `build-mpi-eigenexa-wannier-lib`, shell/Python postprocessing.

---

### Task 1: Capture Current Failure as a Reproducible Diagnostic

**Files:**
- Create: `/private/tmp/c64_owner_local_redef_manifest.tsv`
- Use: `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2c_laser_fraglocal_wfpw_pw128_nt4500_fluxscatterfix`

**Step 1: Record current short-run metrics**

Run existing short variants or create 1000-step variants for:

```text
current global
fragment_local all
fragment_local w_pself
fragment_local w_only
```

Compute:

```text
Pz_ptp_0_100
correlation with A proxy
post-pulse Pz when available
```

**Expected:** `all` has some early correlation but residual/low-frequency shift; `w_only` has almost no response; `w_pself` is weak or phase-wrong.

### Task 2: Add an Owner-Local Field Block Kind

**Files:**
- Modify: `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/src/rt/dg/rt_dg_integrator_expdiag.f90`

**Step 1: Extend input dispatch**

Add a new accepted value for `dg_mixed_z_frag_local_field_block`:

```fortran
case ('owner_local','owner','owned')
  mixed_z_frag_local_field_block_kind = 'owner_local'
```

**Step 2: Map `owner_local` to dynamic block size**

Treat `owner_local` as:

```text
field_block_dim = W_owner_slots + P_self_slots
field_neighbor_slots = 0
```

**Step 3: Verify parsing**

Run a 1-step input using `dg_mixed_z_frag_local_field_block='owner_local'`.

Expected:

```text
field_block_kind=owner_local
field_neighbor_slots=0
replacement_applied=T
bad=F
```

### Task 3: Make Production Backend Use Owner-Local as the Default Fragment Route

**Files:**
- Modify: `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/src/rt/dg/rt_dg_integrator_expdiag.f90`

**Step 1: Ensure owner-local excludes P_neighbor from `nblock`**

In `build_fragment_local_storage_direct`, make `owner_local` use:

```text
p_field_max = 1
```

This keeps `P_owned`/self in the field block and excludes neighbor P from dynamic coefficients.

**Step 2: Preserve diagnostic `all`**

Keep `all` available, but add log wording that it is a diagnostic non-closed route when `P_neighbor_storage_slots > 0`.

**Step 3: Verify no propagated coefficient is dropped**

For `owner_local`, confirm:

```text
P_neighbor_storage_slots > 0 may exist as layout metadata
field_neighbor_slots=0
only P_self contributes to pcoef writeback
```

### Task 4: Short-Run Physics Check

**Files:**
- Create temporary input under `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/`

**Step 1: Run field-off 200 steps**

Expected:

```text
end SALMON
Pz_ptp <= 1e-12
```

**Step 2: Run laser 1000 steps**

Compare owner-local against current global:

```text
Pz_ptp_0_100
correlation with A proxy
Pz trend/DC component
```

Expected:

```text
correlation is not near zero
Pz response is not dominated by monotonic/DC shift
```

### Task 5: Decide Whether Owner-Local Is Sufficient

**Files:**
- Update: `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/docs/plans/2026-07-01-owner-local-mixedz-redefinition-design.md`

If owner-local response is still too weak, document the next step as an embedding correction:

```text
H_eff(owner) = H_owner + boundary correction from neighbor P space
```

The embedding must alter the owner block matrix, not add neighbor P as propagated coefficients.

### Task 6: Commit the Redefinition Checkpoint

**Files:**
- Add only the design/plan docs and focused source changes.

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected:

```text
[100%] Built target salmon
```

Commit message:

```bash
git commit -m "fix: redefine fragment local mixed-z owner block"
```
