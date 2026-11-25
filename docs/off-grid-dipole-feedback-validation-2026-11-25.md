# Off-grid dipole feedback — 2026-11-25

Implemented reciprocal eight-site trilinear macro sampling and dipole-current deposition in Fortran. DipoleMap.from_grid returns the map and zero-based sample indices; MolecularFeedback now accepts multiple sites for one molecule. Existing field-owned Yee sampling/deposition composes with the new map. The own-vacuum reference uses the identical footprint. No numerical integration was added in Python and no extra application-wide communication phase was added.

## Evidence

- Independent affine-field interpolation and the spatial work identity passed at absolute tolerance2e-15. Integrated deposited current equals the molecular dipole increment divided by dt. Grid endpoints, NaN/Infinity and atomic rejection outputs were tested.
- Real vacuum/LD FDTD oscillator tests: outer MPI4,64 steps. Relative midpoint work residuals1.24e-14/1.91e-14; LD scattered feedback3.74348e-9, nonzero. These are molecular work checks, not full-domain energy conservation.
- Existing C2H2 PZ GS reused, no SCF. Position(-.54,1.62,1.94)bohr is off grid. At field27/molecular26, both fields and the molecule were saved; the MPI processes ended. Restoring in new processes and continuing to field64 matched all frames and both FDTD output streams bit-exactly. Restore performed0additional electronic updates. Source and installed-wheel checks passed.
- Previous single-site behavior was compared between the earlier and updated installed wheels: all64frames bit-exact.
- Source related suite63passed. Final sdist-to-wheel installation:92tests,90passed/2optional skips. Manual synchronization and whitespace checks passed; Fortran source/object eigensolver audits passed. No CI launched.

The final wheel SHA256 is `9b14300c1eef548cf8c0acb8246acf878f25eeb6e8a281bc9865a299477e78b9`.
Detailed records and per-file hashes: local OneDrive calc2026/salmonpy-c2h2-nearfield-20261125/offgrid-feedback-validation.json. Cloud synchronization was not independently confirmed. Failed trial directories/logs were retained; one trial placed a source neighborhood into LD and was correctly rejected, then repeated with the complete source neighborhood in vacuum.

## Scope

This is electric dipole feedback with a uniform molecular driving field, not the reciprocal full-current/magnetic return for GridFieldSeries/Iwasa spatial Hamiltonians. That extension, discrete charge-continuity verification and physically converged tip/bow-tie campaigns remain open. Do not interpret these software checks as molecular spectra or a guarantee of physical convergence.
