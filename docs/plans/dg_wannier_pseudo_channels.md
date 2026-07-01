# DG-Wannier Pseudo-Channel Projection Design

## 1. Motivation

The DG-Wannier basis should describe a local low-energy electronic subspace
that remains useful after optical excitation.  It should not be limited to the
occupied valence manifold or to a material-specific projection such as
`C:sp3`.

The target workflow is:

```text
KS orbitals
  -> AO candidate generation
  -> KS subspace selection
  -> Wannier localization
  -> DG diagnostics
  -> export
```

This allows the same machinery to be used for diamond, Si, BN/hBN, water, and
more general materials.

## 2. Current Limitation: C:sp3 Hard-Coded Projection

The current DC-LCFO Wannier export accepts only `wannier_projection='C:sp3'`.
For a C64 diamond cell this creates 256 projection functions and explicitly
stops when `wannier_num_wann < 256`.

That is too restrictive for convergence studies such as:

```text
occ_only
occ + 0.5 Nocc
occ + 1.0 Nocc
```

It also does not transfer to non-carbon systems or systems where polarization
channels are needed.

## 3. New Projection Mode: pseudo_channels

Add a new projection mode:

```text
wannier_projection = 'pseudo_channels'
```

The existing `C:sp3` path remains available as a compatibility and regression
reference.  The new path generates AO-like candidate projectors from the
available pseudopotential channel information and a small default
polarization-shell rule.

The initial implementation should support real Gaussian-like `s`, `p`, and
`d` projectors.  Higher `f` projectors can be added later when needed.

## 4. AO Candidate Generation

For each atom, generate projectors according to pseudopotential angular
channels:

```text
l = 0 -> s
l = 1 -> px, py, pz
l = 2 -> dxy, dyz, dzx, dx2-y2, dz2
```

If the pseudopotential lacks a polarization channel but `include_polarization`
is enabled, add one shell above the maximum available angular momentum, capped
at `d` in the first implementation.

The radial part is initially Gaussian with width controlled by the existing
`wannier_projection_width`.  This keeps the first implementation simple and
close to the existing C:sp3 seed generator.  Later work can use projector norm,
reference energy, or pseudopotential radial tables.

## 5. num_bands / num_wann / num_proj Relation

Use the following relation:

```text
num_wann <= num_bands
num_proj >= num_wann
```

The initial implementation may default to:

```text
num_wann = min(num_proj, num_bands)
```

but it must not require `num_wann >= num_proj`.  In particular, `num_wann <
num_proj` must be valid.

When `num_proj > num_wann`, reduce the projection set by projectability/SVD so
the exported AMN matrix has exactly `num_wann` columns.

## 6. Projectability-Based Pruning

For each candidate projector `p`, compute its norm inside the selected KS
subspace:

```text
q_p = sum_n |<psi_n|p>|^2
```

The first implementation should:

1. Build all AO candidate overlaps.
2. Form the projector overlap/Gram matrix.
3. Use SVD or eigen-decomposition to keep the best-conditioned
   `num_wann` projector combinations.
4. Write AMN using the compressed projector columns.

This avoids material-specific selection rules and permits extra candidate
projectors without over-counting.

## 7. Diagnostics CSV

The existing DG-Wannier+BPW diagnostics are the selection mechanism, not a
hard-coded material rule.  For each candidate Wannier space, emit at least:

```text
case_label
num_bands
num_wann
num_proj
Ncond/Nocc
fsum_avg
C_comm_sum
herm
min_sperp
max_sperp
```

Spectrum diagnostics are run after narrowing the candidate set.

## 8. First Validation: Diamond

Diamond is the first validation case because an existing `C:sp3` result is
available.  The first validation should compare:

```text
C:sp3 existing
pseudo_channels occ_only
pseudo_channels occ + 0.5 Nocc
pseudo_channels occ + 1.0 Nocc
```

Required checks:

```text
num_wann < num_proj does not stop
Z_mix remains Hermitian
fsum_avg improves or remains comparable
C_comm_sum does not diverge
min_sperp/max_sperp remain stable
```

The expected first adoption candidate is the smallest low-energy space whose
diagnostics are saturated relative to the next larger candidate.

## 9. Next Validation: Si, BN/hBN, Water

After diamond:

```text
Si
BN or hBN
water box
```

These validate that the projection generation is not carbon-specific and that
the low-energy Wannier space plus BPW complement behaves as intended across
covalent, polar, and molecular systems.
