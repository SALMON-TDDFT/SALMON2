# Windowed-PW/Wannier Hybrid Design

## Purpose

Build a hybrid DG basis that keeps symmetry-complete, sufficiently localized
Wannier blocks and replaces rejected blocks with real-space-windowed plane waves.
The design must preserve the length-gauge RT formulation without restoring the
old globally replicated WPW machinery.

## Decisions

1. Wannier acceptance and rejection are performed on complete symmetry blocks,
   including conjugate partners. Individual orbitals are never removed alone.
2. Each fragment uses covariant partition-of-unity windows and complete reciprocal
   `G` stars. The GS basis and sparsity graph remain fixed during RT.
3. Windowed-PW candidates are projected out of the retained Wannier space:

   ```text
   p_perp = (I - P_W) p
   P_W = sum_a |w_a><w_a|
   ```

   Retained Wanniers are already orthonormal, so production must not rebuild the
   old global `S_WW^-1 S_WP` coefficient matrix. Projection uses only spatially
   overlapping Wannier rows and bounded owner/neighbor exchange.
4. The projected PW block is not internally orthogonalized. Its sparse metric
   `S_PP` is retained, rank-revealed locally by complete symmetry packets, and
   propagated through the generalized equation

   ```text
   i S dc/dt = [H0 + E(t).Z] c,
   S = diag(I_W, S_PP).
   ```

5. `S`, `H`, and `Z` are assembled from the same hybrid basis and use the same
   WW/WP/PP block graph. In particular, W-P orthogonality does not imply a zero
   W-P position block; `Z_WP` remains part of optical transitions.
6. Full-cell `hpsi` projection is the correctness oracle for static hybrid
   operators. Sparse/window-truncated operators are accepted only after numerical
   comparison with this oracle on small systems.

## Basis Construction

### Wannier block selection

Compute localization and representation receipts for every canonical Wannier
symmetry orbit. Accept a block only when all members meet the localization and
conditioning thresholds. Rejected blocks define the target complement dimension;
the threshold is not allowed to silently change the electron count or break a
degeneracy.

### Windowed-PW packets

Use real partition windows `chi_f(r)` satisfying `sum_f chi_f(r)^2 = 1` on the
periodic grid. For every fragment, generate candidates

```text
p_(f,G)(r) = chi_f(r) exp(i G.r).
```

Add or remove a complete `G` star at a time. Window and `G`-star metadata are
bound to a deterministic fingerprint. Windows are fixed after GS initialization.

### Local Wannier complement

For a candidate owned by fragment `f`, identify retained Wanniers whose bounded
support intersects the window support. Compute their overlaps by row-owned grid
quadrature and subtract the projection. A full-cell diagnostic computes all
overlaps and proves that the omitted tail is below tolerance. No persistent
`N_W x N_P` array is permitted.

### Rank selection

Build sparse `S_PP=<p_perp|p_perp>`. Use packet/block pivoting so an entire
symmetry packet is accepted or rejected. Reject nonfinite, indefinite, or
ill-conditioned packets collectively. Do not form `S_PP^-1/2`; retain the metric.

## Operator Construction

Build WW, WP, and PP blocks for:

- metric `S`;
- field-free Hamiltonian `H0`;
- three position components `Z`;
- velocity/current receipts needed by RT observables.

The same support graph and row/column conventions apply to every operator.
Hamiltonian blocks use the existing full-grid SALMON `hpsi` action, including the
nonlocal pseudopotential, rather than a new hand-derived window Hamiltonian.
Position blocks use one documented periodic/local-origin convention and are
validated for Hermiticity and covariance. Propagation and polarization consume
the identical stored `Z`.

## RT Propagation

The coefficient vector is distributed by fragment-owned Wannier and PW packets.
Each matrix-free application exchanges only graph-neighbor coefficient blocks.
The propagator solves the metric action rather than globally orthogonalizing the
basis. The first implementation reuses a tested generalized Krylov/exponential
action with block-local preconditioning; it does not restore the deleted WPW SCF,
fallback, checkpoint, or production-context layers.

The hybrid basis, metric rank, ownership, and sparsity graph are immutable during
RT. This avoids a time-dependent-basis connection term. Density, current,
polarization, and transition observables include WW, WP, PP, and cross terms.

## Scaling And Memory Contracts

Persistent storage is limited to owned basis rows, local/neighbor sparse blocks,
and bounded communication schedules. Forbidden production allocations include:

- replicated full `N_W x N_P` projection coefficients;
- replicated dense hybrid `N_basis x N_basis` operators;
- all candidate PW grid values at once;
- global all-to-all coefficient exchange per RT step.

With bounded window support, bounded graph degree, and a bounded number of PW
packets per fragment, storage and one operator action are linear in system size.
The implementation reports persistent and transient peak bytes separately, plus
neighbor count and communicated bytes. If the required PW packets grow with
system size, the receipts expose that loss of locality rather than hiding it.

## Failure And Provenance Contracts

Every replicated dimension, tolerance, symmetry catalog, window catalog, and
packet catalog is rank-agreed before shape-dependent collectives. Distributed
row IDs are range-checked and exactly-once owned. Allocation and arithmetic
failures are reduced collectively before return. Fingerprints bind the chain

```text
Wannier selection -> windows/G stars -> complement projection -> S/H/Z
-> RT checkpoint -> observables.
```

Any incomplete symmetry packet, metric rank ambiguity, lost W-P orthogonality,
operator non-Hermiticity, or provenance mismatch fails closed.

## Verification Strategy

1. Unit fixtures verify window partition, complete G stars, local projection,
   sparse metric action, length-gauge orientation, and collective adverse paths.
2. MPI 1/2/4/8 fixtures compare fingerprints and numerical results across uneven
   row ownership and local row permutations.
3. Small-system oracle tests compare local projection with the former dense
   `S_WW^-1 S_WP` formulation and compare sparse WW/WP/PP `S/H/Z` with full-cell
   direct quadrature.
4. Static generalized eigenspectra, degeneracies, symmetry representations, and
   transition matrix elements must agree with the retained full-cell reference.
5. Zero-field RT checks metric norm and energy conservation. Field-driven RT
   checks that propagation and polarization use the same position operator.
6. Si64, an oxide, and water/solution cases report accepted Wannier blocks,
   rejected blocks, PW packets, condition numbers, memory, communication, and
   timing. Neighbor truncation is enabled only after the full-cell oracle passes.

## Explicit Non-Goals For The First Stage

- No adaptive basis changes during RT.
- No sparse density-matrix propagation.
- No restoration of the deleted monolithic WPW production route.
- No claim of linear scaling until measured packet and neighbor counts remain
  bounded in the target material classes.
