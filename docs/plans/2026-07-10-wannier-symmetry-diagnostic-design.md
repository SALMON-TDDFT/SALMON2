# Wannier Symmetry Diagnostic Design

## Goal

Make the Wannier position diagnostic truthful without imposing inversion
symmetry on systems that do not have it.

## Design

The production Wannier import path preserves the selected Wannier gauge.  A
direct AMN gauge may transform or rebuild the position matrix, but it must not
trigger an automatic inversion projection.  The existing
`local_inversion_position` mode remains an explicit Si/C diagnostic option.

The global center permutation reports only whether the current Wannier-center
set is closed under the tested operation.  It is not a complete representation
of a crystal symmetry because Wannier functions may also acquire phases or mix
unitarily.  Operator covariance will ultimately be tested with a representation
matrix reconstructed from the real-space Wannier functions.

Until that representation diagnostic is available, the code must distinguish
an unavailable position residual from a numerical zero.  When the center
permutation exists but no position matrix is supplied, the log reports
`z_odd_available=F` and omits `z_odd_res`.  If the center permutation itself
does not exist, the preceding `perm_ok=F` record terminates that diagnostic.

## Generality

Inversion-specific diagnostics run only when an inversion operation was found
and requested.  A later general diagnostic will apply the same covariance test
to every detected space-group operation, using its rotation matrix for the
three position components.  No symmetry operation is added to the physical
Hamiltonian or observable merely to improve an HHG spectrum.

## Wannier90 site symmetry

The bundled Wannier90 contains the symmetry-adapted Wannier implementation,
but `site_symmetry = .true.` requires a `seedname.dmn` file supplied by the
electronic-structure interface.  Wannier90 does not infer the required band and
Wannier representation matrices from the `.win` file or from an internal
spglib call.  SALMON does not currently write `.dmn`, so enabling the keyword
alone is invalid; the existing probe stops while opening the missing file.

Generating `.dmn` is a separate feature.  It requires the same general
space-group representations needed by the operator-covariance diagnostic and
must be tested independently for the Gamma-only supercell workflow.

The current 576-Wannier C64 route uses deterministic `pseudo_channels`, not the
`C:sp3` random fallback.  The fallback remains relevant only when an sp3
projection supplies fewer functions than `wannier_num_wann`.  Bond-center
projection truncation is likewise a separate basis-selection issue.

## Verification

- A static regression test rejects automatic inversion projection after a
  direct AMN gauge is applied.
- The test requires an explicit unavailable-position log branch.
- The MPI/EigenExa/Wannier build must link successfully.
- An import-only smoke run must show `position sym mode=...` without invoking an
  automatic projection and must not print `z_odd_res=0` for unavailable data.
