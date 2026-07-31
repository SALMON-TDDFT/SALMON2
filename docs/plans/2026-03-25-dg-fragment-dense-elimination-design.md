# DG Fragment Dense Elimination Design

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

**Date:** 2026-03-25

**Goal:** Remove remaining always-resident dense fragment-basis matrices from the DG-RT path wherever block-local fragment data is already sufficient, while keeping correctness and staged validation manageable.

## Scope

This design covers the fragment-block part of DG-RT first:

- `S_mat` / `S_mat_prop`
- `H_mat`
- `momentum_mat`
- SOI complex variants (`S_mat_c`, `S_mat_prop_c`, `H_mat_c`, `momentum_mat_c`)

This design does **not** fully eliminate mixed-basis dense storage in the same phase:

- `H_mat_mixed`
- `S_mat_mixed_prop`
- `S_mat_frag_pw`

Those arrays require a separate representation for fragment-PW cross terms and will be handled in a later phase.

## Current Problems

The codebase already introduced block representations:

- `H_mat_blocks`
- `H_mat_kinetic_blocks`
- `S_mat_blocks`
- `S_mat_prop_blocks`
- `momentum_blocks`

However, dense arrays are still kept alive because several call sites still assume:

- dense `matmul` against `S_mat` or `S_mat_prop`
- dense slicing such as `H_mat(1:n,1:n,ispin)`
- dense overlap propagation broadcasts in `ensure_overlap_prop_available`
- dense temporary rebuilds to feed basis update, unitarity, derivative, and projection paths

As a result, the largest fragment-only matrices remain allocated per rank even when the data model is already naturally block-fragment-local.

## Recommended Approach

Use a staged fragment-block-first conversion.

### Why this approach

- The fragment-only operators already have block storage, so most of the representation work is done.
- It isolates correctness risk from the more complicated mixed-basis representation.
- It reduces resident memory before touching the fragment-PW coupling layer.
- It keeps MPI regression debugging tractable because each phase changes one data family at a time.

## Architecture

### 1. Make block storage authoritative

For fragment-only operators, block storage becomes the canonical runtime state.

Dense arrays should be treated as:

- temporary work buffers only when a dense eigensolver or LAPACK call absolutely requires them
- absent in steady-state propagation, density update, and observables paths

This means new code should read and apply:

- `H_mat_blocks`
- `H_mat_kinetic_blocks`
- `S_mat_blocks`
- `S_mat_prop_blocks`
- `momentum_blocks`

instead of dereferencing dense matrices directly.

### 2. Replace dense linear algebra use sites

The main remaining dense consumers are:

- `rt_dg_integrator_derivative.f90`
- `rt_dg_integrator_unitarity.f90`
- `rt_dg_basis_projection.f90`
- `rt_dg_fragment_basis_update.f90`
- `rt_dg_fragment_basis_update_soi.f90`

These should be migrated to block application helpers or small local dense workspaces assembled only when required.

### 3. Remove dense overlap broadcast as steady-state behavior

`ensure_overlap_prop_available` currently broadcasts dense overlap copies inside `icomm_frag`.

This should be replaced with one of:

- block synchronization from root-owned block state
- explicit root-only dense scratch rebuild followed immediately by block refresh and scratch release

The important property is that subgroup ranks must not need a persistent dense overlap copy.

### 4. Keep mixed-basis dense path for a later phase

`H_mat_mixed` and `S_mat_mixed_prop` remain temporarily because:

- they couple fragment-fragment, fragment-PW, and PW-PW sectors
- several routines still require direct mixed dense access
- a good replacement likely needs a dedicated mixed block/distributed operator API

Deferring this keeps phase 1 realistic and testable.

## Data Flow After Phase 1

Fragment-only propagation and update flow should look like:

1. Hamiltonian and overlap construction produce block data.
2. Operator application paths consume block data directly.
3. Any dense rebuild is local, temporary, and scoped to a routine that truly requires dense LAPACK input.
4. Persistent dense operator arrays are no longer required in the steady-state fragment-only path.

## Validation Strategy

Validation should be staged and conservative.

### Build validation

- full local build after each sub-phase

### Behavioral validation

- compare energy/current/charge against the current branch on a small DG-RT case
- verify no regression for both non-SOI and SOI paths where touched
- run with `nproc_frag > 1` after every MPI-affecting change

### Memory validation

- measure resident memory before and after each phase
- distinguish steady-state allocation reduction from temporary dense scratch usage

## Risks

### Risk 1: hidden dense assumptions

Some routines may look block-ready but still rely on dense symmetry, contiguous slices, or direct `matmul`.

Mitigation:

- convert one consumer family at a time
- add small source-level regression checks for forbidden dense access patterns where practical

### Risk 2: SOI/non-SOI divergence

The non-SOI and SOI modules are structurally parallel but not identical.

Mitigation:

- treat each non-SOI change as requiring an explicit SOI mirror review

### Risk 3: overlap propagation correctness

Changing `S_mat_prop` handling can silently break unitarity stabilization.

Mitigation:

- keep a temporary dense fallback behind a clearly isolated rebuild helper until all consumers are migrated

## Deliverables

Phase 1 should deliver:

- fragment-only `S/H/momentum` block storage as the runtime authority
- dense overlap broadcast removed from the steady-state path
- dense fragment-only matrices no longer persistently allocated unless a routine explicitly opts into temporary dense scratch
- mixed-basis dense storage intentionally unchanged and documented as phase 2
