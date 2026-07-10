# SAWF DMN Design

## Goal

Enable Wannier90 symmetry-adapted Wannier functions for the Gamma-only
DC-LCFO workflow by generating a gauge-consistent `seedname.dmn` file from
SALMON wavefunctions and projection orbitals.

## Scope

The first production route supports the current deterministic
`pseudo_channels` basis with complete real s+p+d shells.  The implementation
must be independent of inversion symmetry and use every valid space-group
operation supplied for the actual structure.  Structures with no nontrivial
symmetry use identity only.  Unsupported or incomplete projection shells fail
before Wannier90 is launched.

The initial implementation targets the existing Gamma-only supercell path:

```text
num_kpts = 1
nkptirr = 1
ik2ir(1) = 1
ir2ik(1) = 1
kptsym(:,1) = 1
```

Multi-k SAWF is outside this change.

## User Interface

The `&dc` namelist gains:

```fortran
wannier_site_symmetry = 'off' ! off, auto, file
wannier_symmetry_tolerance = 1d-6
```

- `off` preserves the existing workflow.
- `auto` requires a SALMON build with spglib and detects operations from the
  lattice, fractional coordinates, and atomic species.
- `file` reads the existing `sym.dat` rotation-plus-translation format without
  enabling density symmetrization.

The generated `.win` contains `site_symmetry = .true.` and
`symmetrize_eps = ...` only after `.dmn` construction and validation succeed.
SAWF is rejected when a frozen inner window is requested because Wannier90 3.1
does not implement that combination.

## Build Configuration

`USE_SPGLIB` is an optional CMake feature.  It uses `find_package` and accepts a
site-provided `spglib_ROOT`; SALMON does not download spglib during the build.
Configuration fails when `USE_SPGLIB=ON` and the library is unavailable.

`wannier_site_symmetry='file'` and `'off'` remain usable without spglib.  This
keeps the Fujitsu/EigenExa build independent of network access and permits an
explicit `sym.dat` workflow on Fugaku.

## Architecture

The SAWF implementation lives in a new module,
`src/gs/dc/lcfo_wannier_sawf.f90`.  `lcfo_flux.f90` supplies the DC geometry,
projection metadata, LCFO state access, and output paths, then calls the SAWF
module.  The new module owns:

- normalized space-group operations in fractional and Cartesian coordinates;
- atom and projection-orbital mappings;
- real s, p, and d rotation representations;
- band-space symmetry matrices reconstructed from real-space LCFO states;
- validation and formatted `.dmn` output.

The existing `src/symmetry/symmetry.f90` parser is exposed through a side-effect
free routine for `file` mode.  SAWF use must not implicitly turn on density or
force symmetrization.

## Symmetry Operations

Each operation stores an integer fractional rotation `W`, fractional
translation `tau`, Cartesian rotation `R`, and an atom permutation.  Every
operation is revalidated against lattice metric, species, and periodic atomic
positions, regardless of whether it came from spglib or `sym.dat`.

Validation requires:

- `W` maps the lattice onto itself and `R` is orthogonal;
- atoms map one-to-one onto atoms of the same species;
- identity is present;
- products and inverses close within the supplied operation set.

Invalid operations cause a fatal diagnostic naming the operation and residual.

## Wannier Representation

For `pseudo_channels`, projection indices are reconstructed in exactly the same
atom and real-harmonic order used to write the `.win` and `.amn` files.

- s transforms as a scalar;
- p transforms with the Cartesian rotation matrix;
- d transforms with the corresponding real l=2 representation;
- the atom permutation moves each complete shell to its symmetry partner.

This produces `D_wann(g)`.  Proper and improper rotations, including inversion,
are handled by the same representation.  Tests compare the analytic p/d
matrices against direct polynomial transformation of real harmonics.

`C:sp3`, `Si:sp3`, and `bond_centers` remain disabled in SAWF mode until their
projection representations have dedicated tests.  No random completion is
allowed in SAWF mode.

## Band Representation

For every operation, SALMON evaluates

```text
D_band(m,n;g) = <psi_m | g psi_n>
```

from the reconstructed Gamma-point LCFO states on the periodic real-space
grid.  Grid points are mapped in fractional coordinates, so translations and
non-Cartesian lattice descriptions do not introduce sawtooth coordinates.
Fragment-local contributions are accumulated and reduced over the DC total
communicator; the implementation must not require a permanent global
wavefunction array.

The overlap matrix is polar-orthonormalized only after reporting its singular
values.  If the selected band subspace is not closed under an operation, SAWF
stops rather than projecting away the leakage.

## Consistency Checks

Before writing `.dmn`, SALMON checks:

- unitarity of every `D_band` and `D_wann`;
- group multiplication residuals for both representations;
- band Hamiltonian covariance;
- projection covariance between `.amn`, `D_band`, and `D_wann`;
- identity-operation residuals;
- complete s+p+d shell counts and ordering.

All residuals and the worst operation are written with a stable
`[DC-LCFO-SAWF]` log prefix.  Any failed check prevents `.dmn` and `.win` SAWF
activation.

## Data Flow

```text
lattice + atoms
    -> auto(spglib) or file(sym.dat)
    -> normalized and validated operations
    -> atom/s+p+d projection mapping -> D_wann
LCFO Gamma states + periodic grid
    -> symmetry-transformed overlaps -> D_band
AMN + eigenvalues
    -> representation consistency checks
validated D_wann + D_band
    -> seedname.dmn
    -> site_symmetry lines in seedname.win
    -> Wannier90
```

## Verification

Unit tests cover identity, inversion, 90-degree rotations, improper rotations,
atom permutations across periodic boundaries, p/d real-harmonic matrices,
`.dmn` ordering, and all fail-fast paths.  A parser round trip reads the
generated file using Wannier90 ordering.

Integration proceeds in three gates:

1. A small synthetic Gamma-only cell passes Wannier90 `site_symmetry` without a
   missing or malformed `.dmn` error.
2. C64 pseudo-channel 576/640 import reports closed representations and produces
   symmetry-adapted centers and position matrices.
3. Field-off stationarity and long-pulse HHG are compared with the current
   direct-global reference; no post-hoc inversion projection is enabled.

SAWF becomes the normal route only after all three gates pass.
