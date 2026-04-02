# DG-SOI Nonlocal PP Design

**Date:** 2026-03-30

**Goal:** Connect the existing spin-orbit nonlocal pseudopotential used by the full SALMON SOI path to the DG-SOI Hamiltonian and time-propagation path, without introducing dense/global fallback routes.

## Context

The full SOI path already prepares and applies spin-orbit nonlocal pseudopotentials:

- `src/atom/pp/prep_pp.f90` builds `uv_so`
- `src/so/update_kvector_so.f90` builds `zekr_uV_so`
- `src/so/pseudo_pt_so.f90` applies the SO nonlocal operator to real-space spinor orbitals

By contrast, the current DG-SOI path only uses the scalar nonlocal projector path:

- `src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90` builds nonlocal matrix elements from `ppg%uV` and `ppg%rinv_uvu`
- no code in `src/rt/dg` currently references `uv_so`, `zekr_uV_so`, or `rinv_uvu_so`

This means the DG-SOI path likely misses the spin-orbit part of the nonlocal pseudopotential, even though the full code path includes it.

## Design Principles

- Match the physics of the existing full SOI path first.
- Stay block-local and fragment-local; do not add dense/global fallback machinery.
- Treat SO nonlocal PP as part of the DG Hamiltonian/operator model, not as an afterthought in observables only.
- Keep the scalar nonlocal and SO nonlocal pieces structurally parallel so later verification is straightforward.

## Recommended Approach

Use the existing full-SOI projector data (`zekr_uV_so`, `rinv_uvu_so`) to build a DG-SOI complex matrix contribution

`H_so_nl = <phi_i | V_so^nl | phi_j>`

directly in fragment basis, then add it to the DG-SOI Hamiltonian blocks and derivative application path.

This is the closest DG analogue of `pseudo_so` and avoids introducing a real-space fallback apply inside DG.

## Data Flow

1. PP initialization continues to build:
   - `ppg%uv_so`
   - `ppg%zekr_uV_so`
   - `ppg%rinv_uvu_so`
2. DG-SOI Hamiltonian construction reads those arrays and evaluates spinor projector overlaps against `phi_frag_c`.
3. The SO nonlocal contribution is accumulated into DG complex Hamiltonian data structures.
4. Time propagation uses the DG Hamiltonian that now already contains SO nonlocal PP.

## Scope

### In scope

- DG-SOI static Hamiltonian construction
- DG-SOI adaptive-basis rebuild path if it reconstructs Hamiltonian operators
- DG-SOI derivative path consistency checks if any cached/operator split requires explicit SO-NL refresh

### Out of scope

- Rewriting the full SALMON SOI nonlocal pseudopotential implementation
- New fallback or dense-operator routes
- SOI current/force validation beyond confirming the Hamiltonian/operator path is connected

## Implementation Shape

### 1. Add a DG-SOI SO nonlocal matrix builder

Add a new routine alongside the existing scalar nonlocal matrix construction in:

- `src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90`

This routine should:

- loop over SO projector channels from `ppg%ia_tbl_so`
- build spinor projector overlaps using `ppg%zekr_uV_so`
- apply `ppg%rinv_uvu_so`
- assemble a complex matrix contribution in the DG basis

### 2. Keep scalar and SO nonlocal PP separate in code, then sum

Do not overload the current scalar `add_nonlocal_pp_matrix` logic. Instead:

- keep the scalar nonlocal path for `uV/rinv_uvu`
- add a separate SO nonlocal path for `zekr_uV_so/rinv_uvu_so`
- sum both into the complex DG Hamiltonian representation

This makes debugging and full-vs-DG comparison easier.

### 3. Prefer block-local accumulation

If the DG-SOI Hamiltonian path is already block-centric, accumulate SO nonlocal contributions into the same block structures used by the rest of DG-SOI. Avoid reintroducing dense-only intermediates unless the existing SOI code path already requires them and there is no local equivalent yet.

## Verification Strategy

### A. Structural checks

- confirm DG-SOI now references `zekr_uV_so` and `rinv_uvu_so`
- confirm scalar-only pseudopotential files still run
- confirm SOI pseudopotential files with `yn_spinorbit='y'` no longer ignore SO-NL in DG

### B. Matrix checks

At `A=0`:

- verify the SO nonlocal contribution is finite and nonzero for a system/pseudopotential known to include SOI
- verify Hermiticity of the full DG-SOI Hamiltonian after adding SO-NL

### C. Physics consistency checks

Compare against the full SOI route for a very small system:

- same pseudopotential
- same `A=0`
- compare low-energy matrix/eigenvalue shifts with and without SO-NL enabled

The goal is not exact equality of every DG/full quantity, but confirming that the DG-SOI path responds to SO nonlocal PP in the same qualitative and near-quantitative direction.

## Risks

### Risk 1: Wrong spinor algebra

`pseudo_so` is not just the scalar projector duplicated for two spins. The DG matrix builder must follow the same spinor coupling structure as the full SOI projector application.

Mitigation:

- mirror `pseudo_so` algebra directly from `src/so/pseudo_pt_so.f90`
- do not infer coupling from scalar code

### Risk 2: k/A dependence mismatch

The full code updates `zekr_uV_so` through `update_kvector_so`.

Mitigation:

- DG-SOI matrix construction must consume the already-updated `ppg%zekr_uV_so`
- do not rebuild its own inconsistent phase convention

### Risk 3: Partial connection

If SO-NL is added only to initial Hamiltonian construction but not to basis-update rebuilds or derivative/operator refresh paths, the run will become internally inconsistent.

Mitigation:

- audit every DG-SOI Hamiltonian rebuild entry point before implementation
- use one shared helper where possible

## Recommendation

Implement a new DG-SOI SO nonlocal matrix builder that mirrors full `pseudo_so` physics using `zekr_uV_so` and `rinv_uvu_so`, keep it separate from the scalar nonlocal builder, and wire it through every DG-SOI Hamiltonian rebuild path.
