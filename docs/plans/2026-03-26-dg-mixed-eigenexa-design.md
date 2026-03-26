# DG Mixed Basis EigenExa Design

## Goal

For `n_pw > 0`, stop carrying fragment-PW coupling explicitly through every RT step.
Instead:

1. assemble the mixed fragment+PW overlap/Hamiltonian only at initialization and basis-update time,
2. orthogonalize the mixed basis with the overlap matrix,
3. diagonalize the resulting standard Hermitian problem with the existing EigenExa path,
4. propagate in the resulting orthonormal mixed basis.

This keeps self-consistent updates exact while moving the expensive mixed-basis algebra out of the inner RT loop.

## Current State

The current mixed path keeps raw fragment basis (`F`) and raw plane waves (`P`) and explicitly stores:

- `S_mat_frag_pw`
- `H_mat_frag_pw`
- mixed overlap solves with `FP/PF`
- mixed Hamiltonian application in `calculate_time_derivative`

This means the RT inner loop still pays for:

- mixed overlap handling,
- mixed Hamiltonian application,
- PW-fragment coupling refresh in stage updates.

The repo already has two useful building blocks:

- mixed matrix assembly in `rt_dg_plane_wave.f90`
- standard dense distributed diagonalization with EigenExa (`eigen_sx`) in DG basis-update and DC-LCFO paths

## Recommended Approach

Use the mixed overlap matrix only to build an orthonormal propagation basis.

### Stage A: Build mixed matrices

At initialization and basis update:

- assemble dense mixed overlap `S_mix`
- assemble dense mixed Hamiltonian `H_mix`
- if needed, assemble mixed momentum blocks later as a separate follow-up

The mixed basis order stays:

- fragment rows/cols first
- PW rows/cols next

### Stage B: Orthogonalize

Convert the generalized EVP

`H_mix c = S_mix c e`

into a standard Hermitian EVP.

Recommended method:

1. diagonalize `S_mix`
2. remove tiny eigenvalues with a threshold
3. build `X = U * diag(s^{-1/2})`
4. form `H_ortho = X^H * H_mix * X`

This is preferred over raw Cholesky because it also gives a clean rank cutoff for near-dependent PW directions.

### Stage C: Diagonalize with EigenExa

Use existing EigenExa standard-problem flow:

- distribute `H_ortho`
- call `eigen_sx`
- recover mixed-basis coefficients with `C_mix = X * V`

This matches the current SALMON/EigenExa integration style and avoids introducing a separate generalized-EVP backend.

### Stage D: Propagate in orthonormal mixed basis

After this transformation, the RT propagation basis is no longer raw `F + P`.
Instead it is an orthonormal mixed basis `X_mix`.

The inner RT loop should then prefer:

- identity overlap in the mixed propagation basis,
- no explicit `S_mat_frag_pw` solve in the hot path,
- Hamiltonian application in the transformed basis.

The initial rollout can keep the old mixed path as fallback while the new path is validated.

## Data Model Changes

Add propagation-basis objects to `s_dg_fragment_rt`.

Recommended additions:

- `mixed_basis_ready`
- `mixed_basis_dim(ispin)`
- `mixed_overlap_evals(:,ispin)` for diagnostics
- `mixed_transform(:,:,ispin)` storing `X`
- `coef_mix(:,:,ispin)` or a clear mapping rule between existing coefficient storage and mixed basis

The raw PW data (`k_pw`, `coef_pw`, `S_mat_frag_pw`, `H_mat_frag_pw`) should remain available for rebuild events, diagnostics, and fallback.

## Execution Flow

### Initialization

1. build fragment basis
2. select PW basis
3. assemble `S_mix`, `H_mix`
4. orthogonalize and diagonalize mixed basis
5. project initial coefficients into orthonormal mixed basis

### Basis Update

1. rebuild fragment basis
2. keep or reselect PW basis
3. reassemble `S_mix`, `H_mix`
4. rebuild orthonormal mixed basis
5. project old coefficients into the new basis

### RT Loop

Use the orthonormal mixed basis by default.
Only fallback to raw mixed `FP/PF` handling if mixed re-diagonalization was not prepared.

## Error Handling

- If `S_mix` has too many tiny eigenvalues, emit a clear diagnostic and drop dependent directions.
- If the retained subspace dimension changes after a basis update, log the change and project coefficients conservatively.
- If EigenExa is unavailable, use the existing LAPACK/ScaLAPACK fallback on the standard-form matrix.

## Testing

1. initialization-only mixed basis:
   - verify `X^H S X ≈ I`
   - verify Hermiticity of `H_ortho`
2. basis-update rebuild:
   - verify coefficient reprojection does not explode norm
3. RT smoke test with `n_pw > 0`:
   - compare electron number and total energy trend against current path
4. performance:
   - confirm reduced time in mixed overlap / mixed refresh sections of stage update

## Recommendation

Implement in phases:

- Phase 1: initialization-time mixed orthogonalization + diagonalization, with runtime fallback
- Phase 2: basis-update-time rebuild and reprojection
- Phase 3: remove hot-path dependence on `S_mat_frag_pw` in the common mixed RT path
