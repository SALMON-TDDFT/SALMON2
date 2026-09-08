# DG Basis Update Local Projection Design

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

**Date:** 2026-03-25

**Goal:** Replace the current basis-update-time global re-diagonalization path with a local projection transfer that uses fragment-fragment and fragment-PW overlaps, so RT basis updates no longer require global dense mixed operators to remain resident.

## Problem

The current basis update path still pulls dense mixed operators back into the code because it treats basis refresh as a global eigenbasis reconstruction step:

- update fragment basis functions
- rebuild mixed `H`/`S`
- solve a dense generalized eigenvalue problem
- project propagated coefficients into the refreshed global eigenbasis

This makes basis update the remaining justification for:

- `H_mat_mixed`
- `S_mat_mixed_prop`
- dense mixed work arrays that scale with `n_mat_max + n_plane_waves`

For the current large cases, that memory cost is now the blocker. We cannot verify the old projection strategy if the run dies before or during the update.

## Design Decision

Use **local projection transfer** as the primary basis-update mechanism.

The update should be treated as:

- old fragment basis + old PW basis
- new fragment basis + new PW basis
- coefficient transfer by overlap

and **not** as:

- rebuild a new global eigenbasis for the whole system at each update

This means the runtime object we need during update is an overlap transfer operator, not a dense global Hamiltonian.

## Required Information

For a basis update, the propagated state transfer should use only:

- `S_ff(new, old)` for fragment-fragment overlap
- `S_fp(new, old_pw)` for fragment-to-PW overlap
- `S_pf`
- `S_pp`

The PW sector already represents delocalized physics. That means a basis-update-time transfer does not need remote fragment eigenvectors from far-away fragments just to preserve long-range character. Long-range amplitude is carried through the PW coefficients.

The resulting design assumption is:

- local fragment refresh is local or near-neighbor in the fragment sector
- nonlocal/delocalized content is handled through the PW sector

## Architecture

### 1. Split basis update into two distinct concepts

Keep the following separate:

- **basis refresh**
  - update fragment basis functions
  - refresh fragment operators
- **state transfer**
  - move propagated coefficients from old basis to new basis

The current code mixes these with a global diagonalization step. The new path removes that coupling.

### 2. Make projection the canonical state-transfer path

For basis update, compute the new coefficients from overlap relations:

- fragment-only path:
  - solve for `C_new` using new-vs-old fragment overlap
- mixed path:
  - solve for `C_new_mixed` using `FF/FP/PP` overlap blocks

This should become the normal RT update transfer path. Any dense generalized EVP becomes optional debugging or legacy fallback, not required runtime state.

### 3. Restrict dense work to local temporary solves only

Dense matrices may still appear in tightly scoped local tasks:

- small `nstate_frag x nstate_frag` projection solves
- temporary dense overlap work assembled only for one routine

But there should be no requirement to keep full mixed dense operators resident across RT steps or basis updates.

### 4. Remove basis-update dependence on `H_mat_mixed`

`H_mat_mixed` should no longer be required for:

- coefficient transfer during basis update
- PW compaction bookkeeping
- mixed basis refresh state handoff

If a routine still uses `H_mat_mixed`, it should be because it is explicitly solving a dense mixed eigenproblem in a bounded fallback path, not because normal RT update requires it.

## File-Level Changes

Primary files:

- `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90`
- `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_basis_projection.f90`
- `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`
- `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`

Expected refactor:

- add overlap-based mixed projection helpers
- replace dense mixed diagonalization as the default update handoff
- make PW compaction rebuild from `FF/FP/PP` instead of cached dense mixed operators

## Risks

### Risk 1: projection conditioning

The new transfer solve may become ill-conditioned if overlap blocks are nearly singular.

Mitigation:

- keep regularization local
- use explicit condition checks
- zero-cost fallback only if the local solve fails

### Risk 2: behavior change versus old global eigenbasis transfer

The new path preserves propagated state by overlap transfer, not by re-expressing the system in a fresh global eigenbasis.

Mitigation:

- treat this as an intentional algorithmic change for RT scalability
- validate charge/norm stability first
- defer apples-to-apples physics comparison until runs are stable again

### Risk 3: mixed update bookkeeping still hidden in old dense arrays

Some parts of PW compaction or update staging may still assume dense mixed caches exist.

Mitigation:

- explicitly remove `H_mat_mixed` / `S_mat_mixed_prop` from normal update contracts
- rebuild any necessary small slices directly from `FF/FP/PP`

## Success Criteria

This phase is successful when:

- normal RT basis updates no longer require persistent `H_mat_mixed`
- normal RT basis updates no longer require persistent `S_mat_mixed_prop`
- basis-update-time state transfer uses overlap-based local projection
- large runs can pass the basis update stage without dense-memory blow-up
