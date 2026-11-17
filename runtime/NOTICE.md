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
