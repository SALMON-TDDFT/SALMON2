# DG-RT Fragment-Local Orbital Parallelization Specification

Date: 2026-05-04
Status: Draft
Owner: DG-RT development

## Background

The current DG-RT fragment-local MPI layout uses one fragment subgroup per
fragment.  Inside each subgroup, ranks are treated as real-space ranks through
`id_frag`, `isize_frag`, and `nproc_rgrid`.  This has repeatedly made density,
current, overlap, and mixed fragment-plane-wave paths fragile because the code
must keep three ownership models consistent at once:

- real-space ownership for fragment buffer/core regions
- basis-row ownership for fragment basis blocks
- orbital/state ownership for density and coefficient contractions

Recent `n_pw > 0` debugging showed that the physically relevant invariants are
coefficient-space invariants such as raw electron count, `C^H S C`, current,
and `Tr[C^H H C]`.  The normalized real-space density can hide coefficient
drift.  Therefore the fragment-local parallelism should be moved to an
orbital/state decomposition, and real-space decomposition inside a fragment
subgroup should be removed from the DG-RT execution path.

## Goals

- Make fragment-local DG-RT parallelism orbital/state based.
- Keep fragment-level MPI decomposition: one fragment subgroup per fragment.
- Make all ranks inside a fragment subgroup see the same fragment real-space
  buffer, basis data, overlap, Hamiltonian, momentum, and mixed PW matrices.
- Split occupied/state columns across ranks inside the fragment subgroup.
- Reduce orbital partial sums for density, current, energy, and diagnostics.
- Preserve coefficient-space invariants without relying on density
  renormalization.
- Make the C64 2x2x2 diamond, 256-electron, `np=8` case the primary
  acceptance target.

## Non-Goals

- Do not keep fragment-local real-space parallelism as a production path.
- Do not optimize memory first.  The first implementation may replicate
  fragment-local matrices inside each fragment subgroup.
- Do not redesign the global fragment decomposition.
- Do not change the physics of the mixed fragment/PW basis.
- Do not use normalized `rho%f` as the primary correctness metric.

## Parallel Layout

The global layout remains:

```text
np = n_frag * nproc_orb
```

where `nproc_orb` is the number of MPI ranks assigned to each fragment subgroup.
For compatibility, the first implementation may continue reading
`product(nproc_rgrid)` as `nproc_orb`, but the meaning changes in DG-RT orbital
mode:

```text
ifrag_group = global_rank / nproc_orb + 1
id_frag     = global_rank mod nproc_orb
isize_frag  = nproc_orb
icomm_frag  = communicator for one fragment subgroup
```

In orbital mode, `id_frag` is not a real-space rank.  It is the orbital/state
rank inside one fragment.

## Real-Space Ownership

All ranks in a fragment subgroup use the same real-space fragment region.

- `rank_core_lo/rank_core_hi` are identical for all ranks in the subgroup.
- `rank_buf_lo/rank_buf_hi` are identical for all ranks in the subgroup.
- Fragment buffer basis arrays are replicated inside the subgroup.
- Real-space subgroup send/recv/alltoall paths are disabled in orbital mode.
- Density ownership is not assigned by real-space rank inside a fragment.

This removes the need for fragment-local real-space stitching in density and
Hamiltonian updates.  Only orbital partial sums are reduced.

## Orbital Ownership

Each fragment subgroup partitions state columns:

```text
nocc_per_rank = ceil(nocc_spin / isize_frag)
io_s_frag     = id_frag * nocc_per_rank + 1
io_e_frag     = min((id_frag + 1) * nocc_per_rank, nocc_spin)
```

The same partitioning rule is used for:

- density reconstruction
- raw charge and `rho_global_raw_elec`
- current
- energy
- `C^H S C` diagnostics
- mixed-basis occmap diagnostics
- occupied subspace reorthonormalization

For routines that operate on all propagated states, use `nstate_tot` instead
of `nocc_spin` with the same block rule.

## Data Ownership

### Replicated Per Fragment Subgroup

The following data are replicated across ranks inside `icomm_frag`:

- fragment basis and buffer basis
- `index_basis`
- `S_mat_blocks`, `S_mat_prop_blocks`
- `H_mat_blocks`, kinetic/local/nonlocal block data
- momentum blocks
- `S_mat_frag_pw`
- `H_mat_frag_pw`
- `P_mat_frag_pw` and G-space variants
- `H_mat_pw`, `H_mat_pw_diag`
- mixed transform and mixed eigenvalues

### Distributed By Orbital/State Column

The following data are logically distributed by state column:

- `coef(:, io, ispin)`
- `coef_pw(:, io, ispin)`
- `coef_mix(:, io, ispin)`
- RK stage arrays `k`, `k_pw`
- temporary wavefunction batches

The first implementation may keep full coefficient arrays allocated on all
ranks and update only owned columns, then synchronize.  After correctness is
established, memory can be reduced by storing only owned columns plus a gather
cache where required.

## Synchronization Rules

Coefficient synchronization is mandatory after any operation that updates
state columns.

Required sync points:

- after initial coefficient read/projection
- after each RK/SSPRK/AETRS intermediate stage update
- after final time-step update
- after mixed-basis projection or reorthonormalization
- before any routine requiring all columns, such as dense mixed projection
  or diagnostics over the full occupied subspace

The synchronization primitive should be a column-wise allgather or allreduce
over `icomm_frag`.  It should not use real-space rank ownership.

## Density Reconstruction

Density reconstruction becomes an orbital reduction:

1. Each rank reconstructs density from its owned occupied states.
2. Each rank accumulates unnormalized local density contribution.
3. `icomm_frag` reduces density contributions inside the fragment subgroup.
4. Fragment-level/global reductions combine fragment contributions.
5. `rho%f` may still be normalized to the target electron count, but this is
   not a coefficient correctness check.

The raw charge before scaling must be retained and reported:

```text
rho_global_raw_elec = integral of unscaled density contribution
```

This value is the primary density-side invariant for debugging coefficient
drift.

## Current And Energy

Current and energy must use the same orbital ownership as density.

For each owned state:

```text
J_rank += occ(io) * <psi_io | J_op | psi_io>
E_rank += occ(io) * <psi_io | H    | psi_io>
N_rank += occ(io) * <psi_io | S    | psi_io>
```

Then reduce inside `icomm_frag` and across fragment groups as needed.

Important requirements:

- Do not compute energy from normalized density when diagnosing propagation.
- `E_one` must be compared against raw coefficient norm and `C^H S C`.
- Diamagnetic current should use the physical electron count, not the
  accidentally scaled raw coefficient count.
- Debug output must keep both raw and scaled charge.

## Mixed Fragment/PW Path

For `n_pw > 0`, the mixed path must obey one metric convention:

- `mixed_transform` columns represent S-orthonormal mixed modes.
- `coef_mix` stores coordinates in that mixed basis.
- raw coefficients are reconstructed by `raw = T * coef_mix`.
- mixed coefficients are reconstructed from raw coefficients by
  `coef_mix = T^H * S * raw`.

No phase or metric operation may be applied to only one side of this relation.
Any phase adjustment of `T` must be compensated in `coef_mix` or the density
matrix.

The mixed path must synchronize columns after:

- `sync_mixed_coef_from_raw`
- `sync_raw_coef_from_mixed`
- mixed occupied-subspace reorthonormalization

## Removed Or Disabled Paths

In orbital mode, the following fragment-local real-space paths are not used:

- subgroup real-space density send/recv maps
- subgroup real-space alltoallv density stitching
- real-space rank-specific fragment buffer ownership
- basis-row ownership that depends on `id_frag` as a real-space rank
- `coef_owner` logic used to represent basis-row ownership inside a fragment

If any of these paths are reached in orbital mode, the code should fail early
with a clear diagnostic rather than silently mixing ownership models.

## Input And Configuration

Add an explicit DG fragment parallel mode:

```fortran
dg_fragment_parallel_mode = 'orbital'
```

Allowed values:

- `'orbital'`: production mode; fragment-local state-column parallelism
- `'legacy_realspace'`: optional temporary compatibility mode during migration

After migration, `'orbital'` should become the default and legacy real-space
mode should be removed or kept only behind a developer/debug flag.

For the first implementation, `product(nproc_rgrid)` may define
`nproc_orb`.  Longer term, add a clearer input name:

```fortran
nproc_orbital_dg = ...
```

## Diagnostics

The following diagnostics are required for acceptance:

- `rho_global_raw_elec`
- `rho_global_scaled_elec`
- `rho_coef_sum_elec`
- `rho_coef_minus_raw_elec`
- `coef_occ_norm`
- `csc_occ_identity_norm`
- `csc_occ_identity_max`
- `E_one`
- `E_tot`
- current components
- `same_ref1` and transition diagnostics for mixed occmap runs

Correctness must be judged from raw coefficient-space values.  Normalized
`rho%f`, `rho2`, and `rhomax` are secondary diagnostics only.

## Migration Plan

### Phase 1: Layout Gate

- Add `dg_fragment_parallel_mode`.
- Introduce orbital mode initialization.
- Keep one fragment subgroup per fragment.
- Make subgroup ranks share identical real-space fragment buffer metadata.
- Fail early if legacy real-space subgroup communication is called in orbital
  mode.

### Phase 2: Orbital Reductions

- Route density reconstruction through state-column ownership.
- Route current and energy through the same ownership.
- Route raw charge diagnostics through the same ownership.
- Verify `n_pw=0` against the current stable buffer-basis result.

### Phase 3: Coefficient Synchronization

- Add column synchronization after each time integrator stage.
- Apply the same synchronization to `coef`, `coef_pw`, and `coef_mix`.
- Verify `C^H S C` and `rho_global_raw_elec` stability.

### Phase 4: Mixed PW Stabilization

- Verify the mixed fragment/PW metric convention.
- Keep phase-fix disabled unless coefficients are gauge-compensated.
- Verify `n_pw > 0` on H2 as a small regression case.

### Phase 5: Production Acceptance

- Run C64 2x2x2 diamond, 256 electrons, `np=8`.
- Compare `n_pw=0` and `n_pw>0`.
- Confirm no drift in raw charge, coefficient metric, current, or energy.
- Make orbital mode the default.

## Acceptance Criteria

The implementation is accepted when all of the following pass:

- Build succeeds.
- `n_pw=0` buffer-basis DG-RT remains stable.
- `n_pw>0` H2 small regression runs without crash.
- H2 `rho_global_raw_elec` no longer shows artificial coefficient blow-up.
- C64 2x2x2 diamond, 256 electrons, `np=8`, buffer-basis run is stable.
- `C^H S C` remains bounded and close to identity for occupied states.
- `E_one` and current do not diverge while normalized density remains stable.
- Results are independent of fragment-local orbital rank count within expected
  numerical tolerance.

## Risks

- Replicating matrices inside each fragment subgroup increases memory.
- Column synchronization after every RK stage can add communication overhead.
- Some existing code paths may assume `id_frag` is a real-space rank.
- Legacy diagnostics may accidentally combine raw and scaled charge.
- SOI/non-collinear paths may need a separate migration pass.

## Open Questions

- Should the short-term input reuse `nproc_rgrid` or immediately introduce
  `nproc_orbital_dg`?
- Should legacy real-space mode remain available during the migration window?
- What tolerance should be used for C64 `C^H S C` and raw charge conservation?
- Should matrix construction be replicated first, or should some matrix-block
  construction remain distributed after orbital mode is stable?

