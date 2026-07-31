# DG Flux Math Audit

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## Current Symptom

For the 2x2x2 diamond DC-LCFO case:

- no-Flux LCFO gap is about 4.57 eV.
- Flux-LCFO block-diagonal gap is about 0.70 eV.
- Projecting the Flux Hamiltonian into the no-Flux LCFO subspace gives:
  - full face-neighbor Flux: about 0.19 eV
  - current block-diagonalized Flux: about 1.35 eV when using no-Flux coefficients

This points to the Flux Hamiltonian construction, before Wannier90/AA_R, as the
primary source of the collapsed gap.

## SIPG Reference Form

For the kinetic operator

```text
T = -1/2 Delta
```

a symmetric interior-penalty DG form on an interface F between left L and right R
with normal n pointing from L to R can be written as

```text
a(u,v) =
  sum_K 1/2 int_K grad(v) . grad(u)
  - 1/2 sum_F int_F {d_n u} [v]
  - 1/2 sum_F int_F {d_n v} [u]
  + eta sum_F int_F [v] [u]
```

with

```text
[u]       = u_L - u_R
{d_n u}  = 1/2 (d_n u_L + d_n u_R)
```

where `d_n u_R` is differentiated along the same interface normal n, not along the
right element's outward normal.

Expanding one interface gives the matrix factors:

```text
LL: -1/4 v_L d_n u_L -1/4 d_n v_L u_L + eta v_L u_L
LR: -1/4 v_L d_n u_R +1/4 d_n v_L u_R - eta v_L u_R
RL:  1/4 v_R d_n u_L -1/4 d_n v_R u_L - eta v_R u_L
RR:  1/4 v_R d_n u_R +1/4 d_n v_R u_R + eta v_R u_R
```

The current signs in `lcfo_flux.f90` and `rt_dg_fragment_hamiltonian.f90` match
the LL/LR expansion if the remote derivative is already converted to the receiver
normal. So the first suspect is not the simple sign pattern.

## Implementation Mismatch Risks

### 1. Strong-form volume term

`lcfo_flux.f90` builds the volume block as:

```fortran
mat_H_local(io,jo) = sum(f_basis(:,:,:,io) * hf(:,:,:,jo)) * hvol
```

where `hf` comes from applying the regular FD `hpsi` stencil to a basis function
that has been written only on the fragment core.

The SIPG formula above assumes either:

- a weak volume term `1/2 int_K grad(v).grad(u)`, or
- a strong term `int_K v (-1/2 Delta u)` without zero-extension boundary artifacts.

Applying the SALMON FD stencil to a core-truncated function can already encode a
Dirichlet-like boundary contribution at fragment faces. Adding SIPG surface terms
on top of that may double count or distort the kinetic boundary operator.

This should be tested by writing separate kinetic components:

```text
H_volume_strong
H_surface_self
H_surface_cross
```

and checking the gap as these are combined.

### 2. Block-diagonalization is not the DG global operator

The full DG Hamiltonian is block-sparse:

```text
H_DG = H_II + H_IJ(face neighbors)
```

The current `SALMON_DG_BLOCK_DIAG_H` path folds face-neighbor information into
fragment-local blocks through symmetrized halo terms. This is not algebraically
equivalent to the full DG operator unless an additional approximation or
downfolding argument is supplied.

In particular, off-diagonal Flux terms multiply neighbor coefficients. Replacing
them with local diagonal-block contributions changes the eigenproblem.

### 3. Interface counting in 2x2x2 periodic cells

For two fragments along an axis, the + and - periodic neighbors can map to the
same fragment. Interface contributions must still be counted as distinct physical
faces, but matrix folding and 0.5 factors must be checked carefully because the
same `(I,J)` block can receive two geometrically different faces.

### 4. Penalty parameter scale

The current value is:

```fortran
alpha = 10 / h
```

This may be numerically large for high-order FD-derived basis functions. A
penalty sweep should be diagnostic only, not a final fix, but it can reveal
whether the gap collapse is dominated by the penalty term.

## Next Diagnostic Steps

1. Add optional component output in `lcfo_flux.f90`:
   - strong volume core
   - SIPG self surface
   - SIPG cross surface
   - velocity/position Flux separately

2. Use `tools/dg_lcfo_flux_gap_projector.f90` to evaluate gaps for:
   - volume only
   - volume + self surface
   - volume + full surface
   - block-folded surface

3. Repeat with a kinetic weak-form volume matrix:

```text
1/2 int_K grad(v).grad(u)
```

instead of `int_K v hpsi(u)` and compare gaps.

4. Only after the static spectrum is physically reasonable should Wannier90 and
   AA_R be brought back into the path.

## 2026-06-18 Component Projection Result

The component dump was added as `hamiltonian_flux_components.bin` and analyzed
with `tools/dg_lcfo_flux_gap_projector.f90`.  The diagnostic fixes the
coefficient subspace to the no-Flux LCFO eigenvectors, so the comparison probes
the Hamiltonian construction rather than a changed basis.

For the 2x2x2 diamond case with `nstate=512`, `nocc=128`:

```text
reference no-Flux H in no-Flux LCFO subspace: 4.568937 eV

Flux-H components projected to no-Flux LCFO subspace:
  volume only       : 1.297733 eV
  surface self only : ~0 eV
  volume + self     : ~0 eV
  cross only        : 8.616640 eV
  full sparse H     : 0.189204 eV
  block-folded H    : 1.348491 eV

Flux-H block eigenvectors with block-folded H:
  gap               : 0.697497 eV
```

This rules out Wannier90, AA_R, length-gauge matrix elements, and RT time
propagation as the primary cause of the collapsed gap.  The gap is already
unphysical in the static LCFO-Flux Hamiltonian.

The most important observation is that `volume only` is already far below the
no-Flux reference.  Therefore the first root-cause target should be the volume
kinetic matrix construction:

```fortran
sum(f_basis(:,:,:,io) * hf(:,:,:,jo)) * hvol
```

where `hf` is obtained by applying the regular SALMON `hpsi` stencil to a
fragment-core-truncated basis.  This is not the same discrete operator as the
SIPG weak form volume term

```text
1/2 int_K grad(v) . grad(u)
```

and can already include implicit zero-extension/Dirichlet-like face artifacts.
Adding explicit DG surface terms on top then combines incompatible boundary
operators.

### Immediate Conclusion

The Flux sign pattern should not be changed first.  The next implementation
step should be a diagnostic weak-form kinetic volume matrix, using the already
available local basis gradients, and then the same projection table should be
recomputed:

```text
H_volume_weak = 1/2 sum_axis int_K grad_axis(v) grad_axis(u)
             + local/pseudo/Hartree/XC potential terms
```

If the weak-form volume restores the no-Flux-scale gap before adding surface
terms, the production Flux Hamiltonian should be rebuilt around that weak
volume form.  The block-folded Hamiltonian remains a separate approximation and
should only be used after the sparse DG operator itself has a physical static
spectrum.

## 2026-06-18 Weak-Volume Diagnostic

A second diagnostic file, `hamiltonian_flux_weak_components.bin`, was added.
It keeps the present local/pseudo/Hartree/XC contribution from `hpsi`, but
replaces the core-zero kinetic part by a buffered weak-form kinetic estimate:

```text
H_volume_weak = H_volume_strong
              - T_core_zero_strong
              + 1/2 int_K grad(v_buffered) . grad(u_buffered)
```

Projection to the same no-Flux LCFO subspace gives:

```text
weak volume only       : 0.435794 eV
weak volume + self     : ~0 eV
weak full sparse H     : ~0 eV
weak block-folded H    : 0.114835 eV
```

This does **not** restore the physical gap.  Therefore the problem is broader
than a simple replacement of the kinetic volume by a weak gradient integral.
The surface self term is especially dangerous: both strong and weak volume plus
self surface become nearly degenerate at the HOMO/LUMO boundary.

### Updated Root-Cause Assessment

The simple LL/LR sign pattern remains internally consistent, but the actual
operator being diagonalized is not a reliable DG discretization of the original
Hamiltonian.  The likely issues are:

1. the LCFO basis is first localized/truncated/orthogonalized on fragment cores,
   so the local block spectrum no longer resembles the no-Flux LCFO spectrum;
2. the self surface term is being added as a large local boundary operator
   without a corresponding variationally consistent volume/boundary pairing;
3. the block-folded approximation is not algebraically equivalent to the sparse
   DG Hamiltonian and further changes the eigenproblem.

### Practical Next Step

Before using this Flux Hamiltonian for Wannier or RT, build a static
verification path that checks the following in order:

```text
1. fragment-local volume-only LCFO spectrum versus no-Flux LCFO reference
2. volume + self surface spectrum
3. full sparse face-neighbor DG spectrum
4. only then any block-folded/downfolded local approximation
```

The current implementation fails already at steps 1 and 2.  The production path
should therefore keep the no-Flux LCFO/DC Hamiltonian for physical validation
and treat the current Flux-H path as experimental until the static spectrum is
repaired.
