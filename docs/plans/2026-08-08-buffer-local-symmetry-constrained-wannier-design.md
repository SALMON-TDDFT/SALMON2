# Buffer-Local Symmetry-Constrained Wannier Design

## Objective

Generate localized orbitals without materializing all full-system Kohn--Sham
wavefunctions or a dense full-system Wannier gauge.  The result must preserve
the exact symmetry of the instantaneous atomic structure, minimize a periodic
spread functional inside overlapping buffered fragments, and converge as the
buffer is enlarged.

The result is called a *buffer-converged symmetry-constrained localized
Wannier basis*.  It is not called an MLWF at finite buffer.  Wannier90 is an
optional small-system oracle, not a production dependency of this route.

## Architecture

The existing complete-s+p projection and rank-revealing construction supplies
only the initial gauge.  A new localization stage applies small unitary Jacobi
rotations between Wannier functions whose buffered supports overlap.  It never
forms a dense full-system wavefunction array.

For every exact full-system symmetry operation `g`, a sparse anti-Hermitian
pair seed `K` is promoted to the symmetry-compatible generator

```text
K_sym = |G|^-1 sum_g D(g) K D(g)^dagger
U_sym = exp(K_sym).
```

This group average is required because `D(g)` can mix a complete local orbital
multiplet and need not be a Wannier permutation.  The exponential is evaluated
only on the connected support block generated from the buffered pair graph;
no dense full-system gauge is formed.  The optimizer therefore moves only
inside symmetry-compatible gauges.  If the
instantaneous structure has only identity symmetry, every fragment orbit is a
singleton and no parent-crystal symmetry is restored.

## Data Flow

1. Construct complete-s+p projected local candidates on each buffered
   fragment.
2. Form the retained overlapping basis and its exact fragment-orbit maps.
3. Evaluate periodic localization matrices from the existing phase links,
   weights, and buffered values.
4. Enumerate overlapping Wannier pairs and their connected support blocks.
5. Optimize a sparse pair seed, group-average its anti-Hermitian generator with
   the exact dense Wannier representations, and exponentiate the resulting
   connected block.
6. Accept an iteration only when the periodic spread does not increase and the
   metric orthogonality and symmetry closure gates remain satisfied.
7. Publish the optimized basis, then assemble and project `S`, `H`, `X`, and
   `V` as before.

## Localization Functional

Use the periodic position phases already available on the real-space grid.
For normalized Wannier function `n`, define

```text
z(n,a) = <w_n | exp(i G_a r_a) | w_n>
Omega = sum(n,a) weight(a) * (1 - |z(n,a)|^2).
```

This functional is bounded, origin independent on the periodic cell, and can
be evaluated from buffered tails.  Pairwise unitary rotations are accepted by
an explicit line search.  The implementation must report both `Omega` and the
maximum accepted pair-gradient magnitude.

## Numerical and Physical Gates

- Spread is finite, nonnegative, and monotonically non-increasing.
- The final pair-gradient is below the configured tolerance.
- `W^dagger S W = I` is preserved within the metric tolerance.
- Symmetry representation unitarity and group closure are unchanged.
- Every accepted generator is anti-Hermitian and commutes with each retained
  exact representation within tolerance.
- Boundary value and gradient gates remain satisfied.
- Increasing the buffer produces converged spread, operator matrices, and
  polarization-derived spectra.
- In an exactly inversion-symmetric reference, even-order HHG is suppressed;
  no displaced or thermally disordered calculation is used for that gate.

## Failure Handling

Reject publication rather than silently retaining a non-converged gauge when
the spread rises, the line search stalls above tolerance, a symmetry pair orbit
is incomplete, or orthogonality/closure is lost.  Identity-only symmetry is a
valid case, not an error.

## Scope

The first production implementation retains the accepted regular DC contract:
one buffered fragment per rank and a common retained Wannier rank.  Different
fragment symmetry orbits are supported.  Nonuniform buffer shapes and
different Wannier ranks require variable-size distributed matrix and V3
checkpoint blocks and are a separate extension.

## Wannier90 Oracle

For small fixtures, optionally compare the converged periodic spread and
subspace projectors with a symmetry-adapted Wannier90 calculation using
`site_symmetry = true`.  Production correctness must not depend on
`USE_WANNIER90`.
