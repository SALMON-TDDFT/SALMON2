# Small-gradient HSE formula provenance

The small-gradient WPBEh algebra in the complete replacement
`overrides/hse-explicit-routing/src/xc/hybrid_semilocal.f90` derives from Libxc
5.2.3 `maple/gga_exc/gga_x_wpbeh.mpl` (Copyright (C) 2017 M.A.L. Marques).
The derived formula and constants retain the Mozilla Public License 2.0;
see `licenses/Libxc-MPL-2.0.txt`. SALMONpy's modification combines the E1/log
cancellation analytically in the squared reduced gradient before differentiating.
The surrounding original SALMON/SALMONpy wrapper retains its existing license.
The complete modified source is distributed with SALMONpy; it is not a modification
to the installed Libxc library or the pinned original solver snapshots.

Primary source: https://gitlab.com/libxc/libxc/-/blob/5.2.3/maple/gga_exc/gga_x_wpbeh.mpl
The original Libxc HSE06 mixing convention and screening parameters are retained.

# General-cell scalar FFT-orbit reuse

Complete replacements for `src/xc/exx_k_orbits.f90` and
`src/symmetry/symmetry.f90` retain the bundled fast-PR solver origins/licenses.
The scalar orbit route now uses full-cell gauges, Cartesian k rotations and
wrap phases for general cells, following the existing spinor route. Positive
diagonal cells retain their prior arithmetic. Host-attached zero-impulse response
initialization no longer assumes a z field; actual vectors are checked per step.
The A/B symmetry comparison uses the existing initialization tolerance to avoid
roundoff-only rejections. Original snapshots remain intact; manifests record
complete-file before/after hashes. No installed binary is patched.
