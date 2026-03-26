# DG Mixed Basis EigenExa Plan

## Scope

Introduce an orthonormal mixed fragment+PW propagation basis for `n_pw > 0`, using existing dense mixed assembly and EigenExa standard-problem diagonalization.

## Phase 1: Initialization Path

### 1. Add mixed-basis state to `s_dg_fragment_rt`

- readiness flag
- retained mixed dimension per spin
- transform matrix `X`
- optional mixed-space coefficient storage

### 2. Build standard-form mixed diagonalization routine

New helper flow:

1. assemble dense `S_mix`
2. assemble dense `H_mix`
3. diagonalize `S_mix`
4. build `X = U * diag(s^{-1/2})`
5. form `H_ortho = X^H H_mix X`
6. diagonalize `H_ortho` with EigenExa or existing fallback
7. recover mixed eigenvectors in original basis if needed

### 3. Wire initialization to use the new routine

- run after PW basis selection and mixed matrix assembly
- keep old mixed path available behind a fallback condition

## Phase 2: Basis Update Path

### 4. Rebuild mixed basis after adaptive fragment updates

- refresh `S_mix/H_mix`
- rerun orthogonalization + diagonalization
- invalidate old hot-path mixed caches if basis changed

### 5. Reproject coefficients

- map old coefficient representation into the new orthonormal mixed basis
- verify norm conservation and electron count stability

## Phase 3: RT Hot Path Simplification

### 6. Prefer orthonormal mixed propagation path

In derivative / overlap / stage update:

- bypass raw `S_mat_frag_pw`-driven overlap solves when orthonormal mixed basis is active
- shrink mixed refresh work to rebuild events

### 7. Keep fallback path

Retain current `FP/PF` path for:

- missing EigenExa build
- numerical rank issues
- debugging comparisons

## Verification

- build: `make -C build -j4`
- orthogonality check: `X^H S X`
- compare RT observables for a small `n_pw > 0` case
- time density/stage-update sections before and after

## Risks

- coefficient reprojection mistakes during basis update
- rank truncation changing mixed dimension unexpectedly
- memory growth from storing transforms

## Mitigations

- start with initialization-only enablement
- add explicit dimension and orthogonality diagnostics
- gate basis-update rebuild behind a conservative flag until validated
