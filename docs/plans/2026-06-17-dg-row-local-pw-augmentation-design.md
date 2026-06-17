# DG Row-Local PW Augmentation Design

## Goal

Add plane-wave helper channels to DG-RT without storing full-system coefficient vectors or full fragment-PW operator rows on every rank.

The motivation is the diamond impulse response. The current DC-local basis appears to miss enough high-energy paramagnetic response that `J_tot,z` can remain biased near the diamagnetic current. PW helper channels are still a good direction, but they must not reintroduce global dense state or operator storage.

## Decision

Use fragment-local / row-owner-local PW augmentation.

Each MPI rank may hold:

- coefficient rows it owns,
- fragment blocks for local output rows,
- neighbor or requested coefficient rows needed by local sparse blocks,
- PW helper data only for local or explicitly requested PW rows.

Each MPI rank must not hold:

- an `n_mat_max + n_pw` mixed coefficient vector,
- replicated `coef_pw_full_cache`,
- full `n_mat_max x n_pw` mixed operator blocks for propagation,
- all fragment coefficient rows as a normal RT path.

The current full mixed PW propagation, derivative, density, and unitarity paths are disabled until this design is implemented.

## Architecture

The active DG state remains coefficient based:

```text
psi(t) = sum_i C_i(t) phi_i + sum_p C_p(t) chi_p
```

The implementation stores this state in row-owner form:

- fragment rows use the existing `coef_owner` and `local_coef_global_ids`;
- PW rows use `coef_pw_owner`;
- mixed operators are applied only to output rows owned by the rank.

The row-local operator application is split into four contributions:

```text
Y_F(local) = H_FF(local, needed F) C_F(needed) + H_FP(local, needed P) C_P(needed)
Y_P(local) = H_PF(local P, needed F) C_F(needed) + H_PP(local P, needed P) C_P(needed)
```

The Hermitian relation may be used to avoid storing both `H_FP` and `H_PF`, but the application must still be row-local. If a rank owns PW output rows, it accumulates only those PW rows. Fragment output ranks do not build the full PW output vector.

## Data Layout

The current global-sized arrays are transitional and should not be used by active RT when PW is enabled:

- `S_mat_frag_pw(n_mat_max,n_pw,nspin)`
- `H_mat_frag_pw(n_mat_max,n_pw,nspin)`
- `P_mat_frag_pw(3,n_mat_max,n_pw,nspin)`
- `coef_pw_full_cache(n_pw,nstate,nspin)`

Replace or complement them with local forms:

```text
S_mat_frag_pw_local(n_local_frag_rows,n_pw_needed,nspin)
H_mat_frag_pw_local(n_local_frag_rows,n_pw_needed,nspin)
P_mat_frag_pw_local(3,n_local_frag_rows,n_pw_needed,nspin)
```

For PW-owned output rows, either:

- construct local `H_PF`/`P_PF` rows directly, or
- gather only fragment rows that couple to the owned PW rows and accumulate using conjugate local `H_FP`.

The key rule is that the requested row lists must be derived from sparse local block maps, not from `1:n_mat_max`.

## Density Reconstruction

The current PW density path replicates all PW coefficients through `coef_pw_full_cache`; this is disabled.

The row-local density path should compute real-space wavefunction batches from:

- local fragment coefficient rows already used by the optimized density path,
- only the PW coefficient rows needed for the grid block or fragment handled by the rank.

If the PW helper basis is global plane waves, the density contribution may require a reduction over PW owners. That reduction must be over density/grid blocks, not by first replicating all PW coefficients on every rank.

## Time Integration

Taylor propagation must keep a fixed Hamiltonian during each Taylor operation. Therefore:

- build or update the local mixed Hamiltonian before the Taylor expansion,
- apply the same local mixed operator for every Taylor order,
- do not update Flux or mixed blocks inside the Taylor order loop.

The first implementation target should be Taylor4-PC because it is the current focus. RK paths can follow after the row-local apply is stable.

## Validation

Validation should be incremental:

1. Build-only: PW enabled should stop until the row-local apply is wired.
2. Unit test or smoke case: local-row PW data structures allocate with dimensions tied to owned rows, not `n_mat_max`.
3. Small RT smoke: 2x2x2 diamond with tiny PW count runs without all-row gather stops.
4. Weak-scaling check: compare 2x2x2 and 8x8x8 with the same per-fragment basis and PW helper policy.
5. Physics check: `J_para,z` should cancel `J_dia,z` more strongly; `J_tot,z` should move toward zero-centered damped oscillation.

## Risks

- Truly global plane waves are naturally delocalized. To keep weak scaling, either PW rows must be distributed carefully or the helper basis must become fragment-compatible/localized.
- Nonlocal PP mixed blocks are the highest-risk operator component.
- Density reconstruction can become communication-heavy if PW coefficients are fetched per grid block without batching.
- Overlap conditioning can degrade when PW helpers are close to existing DC basis functions.

## Accepted Direction

Proceed with the row-local/fragment-local PW augmentation path. Do not revive the real-space RT path. Do not use replicated full mixed coefficient vectors as a fallback.
