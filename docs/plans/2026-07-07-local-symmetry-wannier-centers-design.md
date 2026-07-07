# Local Symmetry Wannier Centers Design

## Goal

Ensure DG local and buffer-periodic Wannier seeds preserve the local structural symmetry of the fragment, so centrosymmetric systems do not acquire spurious even-order HHG from an asymmetric effective basis.

## Background

The C64 diamond input is centrosymmetric about `(0.84, 0.84, 0.84)` in the periodic cell. Diagnostics showed that the exported `local_wannier_basis.bin` centers and `buffer_periodic_wannier_basis.bin` centers are not closed under this inversion. The current `bond_centers` option creates physically sensible seed centers, but the subsequent `B^T B` diagonalization can choose arbitrary linear combinations in nearly degenerate subspaces. The resulting orthonormal functions no longer keep the bond-center orbit symmetry.

## Recommended Approach

Build symmetry-aware local Wannier seeds from local structural orbits:

- Detect local equivalent bond-center orbits from the atomic geometry and periodic nearest-neighbor bonds.
- Keep and orthogonalize full orbits rather than individual centers whenever rank pruning is needed.
- Use a projected-position gauge inside the selected subspace, so the final exported centers remain tied to the local bond-center orbit instead of arbitrary overlap eigenvectors.
- Emit diagnostics comparing seed bond-center orbits and final exported centers under the detected local symmetry operations.

This is more general than hard-coding C or Si. For diamond it recovers the four tetrahedral bond centers around each atom and their fragment-local symmetry orbits. For lower-symmetry structures it falls back to smaller orbits, including singletons when no local symmetry is present.

## Alternatives Considered

1. Post-process exported centers by snapping them back to bond centers.
   This is fast, but it can make centers inconsistent with the exported coefficients and position matrix.

2. Detect full crystal space group and symmetrize globally.
   This is clean for perfect crystals, but it does not match the DG fragment-local weak-scaling target.

3. Use local structural orbits and subspace gauge fixing.
   This is the preferred route because it preserves local symmetry while remaining fragment-local.

## Initial Implementation Scope

The first implementation should be diagnostic and minimally invasive:

- Add a DC export diagnostic that prints local bond-center orbit closure and final `wcenter` closure.
- Add a small helper to match centers under periodic inversion and nearest local symmetry operations.
- Do not change production Wannier coefficients until the diagnostic identifies the exact broken layer.

After that, add a gated experimental path for symmetry-gauged local bond-center Wannier generation.

