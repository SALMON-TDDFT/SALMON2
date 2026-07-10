# SAWF Fragment Alignment Design

## Goal

Preserve crystal symmetry in a compressed fragment-local LCFO space by making
the periodic structure, real-space grid, and regular fragment partition share
the same symmetry action before the DC calculation starts.

## Motivation

The C64 `32^3` mesh is divided into eight `16^3` fragment cores.  In the
original coordinates, inversion is

```text
W = -I
tau = (1/8, 1/8, 1/8)
i' = 4 - i  (mod 32)
```

so one source core is split over all eight target cores.  The resulting
fragment-local LCFO direct sum is not inversion closed: the measured aggregate
local leakage is `0.610116`, and the full 664-state LCFO representation has
`closure_residual = 0.99375`.

Iteratively adding all split image directions is not practical.  At a residual
threshold of `1e-3`, the local dimension grows

```text
83 -> 315 -> 483 -> 586 -> 722 -> 890 -> 1047 -> 1173 -> ...
```

and remains unclosed.  This route is therefore rejected.

## Alignment

A uniform translation of every ion in a periodic cell is physically equivalent
in the continuum.  Translating fractional positions by `a` changes a symmetry
translation according to

```text
tau_aligned = tau + (I - W) a  (mod 1).
```

For C64, choose

```text
a = (11/64, 11/64, 11/64)
tau_aligned = (15/32, 15/32, 15/32).
```

The grid action then becomes

```text
i' = 15 - i  (mod 32),
```

so the inversion center is at `7.5` grid spacings inside the first core, not on
a grid point.  More generally, an even-sized core requires its inversion center
to lie half a grid spacing away from the grid-point lattice: `2 c / dx` is an
odd integer.  The resulting index permutation maps each `16^3` core onto one
complete `16^3` core.  It also maps the corresponding buffer box without
cutting it into independent fragment pieces.

## Safety Policy

SALMON must not silently translate user atoms.  Alignment is an explicit input
preparation step that writes a new atom file and a matching symmetry file.  The
original files remain unchanged, and the production input names the aligned
files explicitly.

The alignment tool accepts the cell, mesh, regular fragment counts, atom file,
and normalized symmetry operations.  It searches periodic translations and
accepts one only when all of the following hold:

- every transformed atom maps one-to-one onto an atom of the same species;
- every symmetry operation induces an exact integer permutation of grid
  indices, while symmetry centers may lie at half-grid positions;
- each source fragment core maps wholly onto exactly one target core;
- the mapped buffer geometry is compatible with periodic indexing;
- identity, inverse, and group closure remain valid after translation.

If no compatible translation exists, the tool stops.  It does not fall back to
post-hoc operator symmetrization or polar-unitary repair.

## Interface

The SALMON `&dc` input gains an explicit symmetry-file path:

```fortran
wannier_symmetry_file = 'sym.dat'
```

`wannier_site_symmetry='file'` reads that path.  Existing inputs keep the
historical `sym.dat` default.

The preparation tool is

```text
tools/align_periodic_structure_to_fragments.py
```

and writes atom and symmetry files only when all checks pass.  The C64 formal
sample uses dedicated `atom_sawf_aligned.dat` and
`sym_sawf_aligned.dat` files.

## Runtime Validation

Before constructing `D_band`, SALMON repeats the fragment-map check with the
actual mesh and decomposition.  It logs one bounded line per operation:

```text
[DC-LCFO-SAWF-ALIGN] operation=... fragment_map_ok=T max_targets_per_source=1
```

SAWF stops before seed publication when an operation splits a source core.
The existing rank-local leakage diagnostic remains available as a second,
wavefunction-level check.

## Verification Gates

1. Synthetic tests cover identity, inversion about a half-grid center,
   translated atoms, incompatible partitions, mixed species, and
   multi-operation closure.  They explicitly reject an implementation that
   requires an even-grid inversion center to coincide with a grid point.
2. The aligned C64 atom and symmetry files map all atoms and all grid points
   exactly, with one target fragment per source fragment.
3. C64 DC-SCF converges to the established 86-iteration result.
4. Identity and inversion `D_band` matrices are unitary within the declared
   tolerance without polar repair.
5. Wannier90 reads the generated `.dmn` and completes SAWF localization.
6. Field-off DG propagation is stationary, then the long `dt=2` pulse is used
   to test whether even-order HHG peaks are suppressed.

## Generality

The algorithm depends on the actual symmetry operations and fragment geometry,
not on C, Si, or inversion specifically.  Non-centrosymmetric and locally
asymmetric structures retain their physical operation set.  A structure for
which no periodic translation makes the chosen regular fragment partition
symmetry compatible is reported as incompatible; supporting it requires a
different partition rather than forced symmetrization.
