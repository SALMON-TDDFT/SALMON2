# Direct FDTD field input validation — 2026-11-25

Implemented one-way molecular input from saved Cartesian E/H fields. Fortran performs component-aware trilinear interpolation and the electric path integral of Takenaka, Taketsugu and Iwasa, JCP 154, 024104 (2021), DOI 10.1063/5.0034933. Orbital magnetic coupling is an additional implementation; that paper omits magnetic interactions.

The FDTD molecule position is independent of the molecular integration center. Axes share their orientation; this release supports translation, not rotation. Paths are split at grid-cell faces. No extrapolation, implicit periodic wrapping, time resampling or live feedback is performed. Input histories must cover the electronic grid, stencil halo and nonlocal projector paths.

## Accuracy

Independent Fortran checks on staggered component grids, with 97 off-grid queries, gave maximum absolute E/H/potential error 2.0816681711721685e-17 for affine fields. Magnetic phase, cell-face crossings, NaN/Infinity rejection and unchanged outputs on invalid queries passed.

For a smooth source-free electrostatic test, spacing 0.5/0.25/0.125 bohr gave relative RMS E errors 6.602/1.455/0.3285 percent and potential errors 3.418/0.6464/0.09014 percent. These are interpolation tests, not convergence evidence for physical tips or bow ties.

## Real field and restart

A two-rank 64-cubed LD/FDTD calculation supplied 4624 local sample sites, with explicit half-step H centering. C2H2 was placed at (21.2, 0.7, 0.9) bohr, off the source grid. Existing PZ ground-state data were reused; no new SCF calculation was performed.

Twenty molecular steps were compared with saving at step 7, terminating the process, and restoring into a new process. All nine observables and output streams were bit-exact. Restoration performed zero extra electronic updates. Independent states remained unchanged. Electron loss was zero in this short software trial.

The source-field window starts at 0.303 atomic time; this provenance is retained. Starting a field-free GS at that window is not equivalent to evolving through its omitted illumination history. Physical field, molecular-box and time-step convergence remain separate checks.

## Distribution checks

A Git-index-only source export produced an sdist and wheel. In a fresh environment importing only the installed wheel, 80 related tests ran: 78 passed, 2 optional skips. The independently linked Fortran analytic test and the real C2H2 direct restart also passed against this wheel. Source inventories passed all 8 tests; Japanese/English guide contents and AI JSON are synchronized. Existing fitted-field results (full, electric-only and uniform controls) were bit-exact against the previously validated executable. No CI was launched.

Executable SHA256: `176e9867a253f2e044ff35e95bebdad751f47e12f338b2b40a39fd7b97535b6f`.
Wheel SHA256: `ab75c36fb7e82f57033b2ea1c78667b2e0005f9a9a0c7ee8c165b6ceb0a83b94`.
Source inventory SHA256: `534dee700504c88acdb29cb13832b99092763c859d11cb5b6100991183c8dd12`.
Field NPZ SHA256: `5da4ec9521dd35b82e3325510055298958d599a70754707f698785f09ead7963`.

Detailed records are in the local OneDrive `calc2026/salmonpy-c2h2-nearfield-20261125/direct-grid-final-wheel-molecule-validation` directory. Cloud synchronization was not independently confirmed.

Remaining work: live field streaming, rotation/nonuniform grids and time resampling, converged tip/bow-tie studies, and two-way feedback.
