# Electric path-response operator — 2026-11-25

Implemented the electric path-integral transpose in Fortran as a stateless host calculation component. Signed charge densities at molecular sample points, quadrature volumes, off-grid molecule placement and staggered E-component origins produce polarization/current grids. Returned currents belong directly to the E-component grids; they must not be passed through the macro-site or dipole deposition map again.

Independent Fortran work absolute/relative errors were6.6174449004242214e-24/1.8435588783747862e-16; integrated dipole absolute error1.6263032587282567e-19. Tests cover nonuniform fields, cell-crossing paths, staggered origins and displaced molecules. Domain, NaN, Infinity and arithmetic-overflow failures retain the entire caller output. Existing forward grid-field source is byte-identical to HEAD; the new response code is separate and does not alter TDDFT/FDTD time evolution.

Source selected tests64passed. A Git-index-only export passed source-inventory preparation, built an sdist and wheel, and passed101installed-wheel tests(99passed,2optional skips). Independent Fortran validation linked against the installed wheel also passed. The saved-array preparation example ran successfully from that environment.

Final wheel SHA256: `20a12ae5987b9c25d8bdb5fa3c0d1be64687bdf2c1b7a6900148fa1803b61c64`.
Detailed logs, synthetic input/output and per-file hashes: local OneDrive calc2026/salmonpy-c2h2-nearfield-20261125/electric-path-response/validation.json. Cloud synchronization was not confirmed. No GS calculation or CI was run.

This is not an end-to-end spatial molecular feedback result. Live density transfer, direct FDTD component-current injection, own-vacuum reference integration, magnetic response and discrete charge-continuity validation remain open. No physical convergence or real-material response is claimed.
