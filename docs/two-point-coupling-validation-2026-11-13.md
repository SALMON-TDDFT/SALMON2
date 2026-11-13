# Two-point native FDTD–TDDFT coupling validation, 2026-11-13

Two independent SALMONpy TDDFT workers were connected in a closed loop to SALMON's
unchanged native 1D Weyl-gauge FDTD routines (`weyl_init`, `weyl_calc`). A validation-only
Fortran process bridge calls those routines and the original incident-field routine;
Maxwell propagation was not reimplemented in Python. The reference is ordinary
`multi_scale_maxwell_tddft` with two MPI ranks, one material macro point per rank.
No production numerical routines or public coupling API were changed by this test.

## Configuration and boundary contract

- Apple Silicon Mac, GNU Fortran/OpenMPI, OpenBLAS, Libxc 5.2.3; one OpenMP and
  BLAS thread per process. The serial TDDFT workers have independent state.
- Two adjacent TDDFT macro points at native indices 1 and 2; FDTD indices -4..6,
  spacing 10 a.u., native absorbing boundaries at both x ends, identity rotation.
- 256 steps, dt=.05 a.u. Incident native Acos2: E amplitude .0007 a.u., frequency
  .7 a.u., duration 8 a.u., z polarization. All three A/current components are exchanged.
- Each material point uses the same converged 8-atom Si initial state, 12³ grid,
  2³ k points, 32 states/electrons. PZ and TB-mBJ were tested separately. TB-mBJ
  uses the reconverged state documented in `tbmbj-validation-2026-11-13.md`.
- Native initialization seeds the incoming pulse and advances the electromagnetic
  field to step zero. EM current starts at zero, as in main_ms.
- Step n: FDTD uses charge current from step n-1; then TDDFT takes A/c at n and
  returns electronic matter current J(n). Store charge current = -J(n) for the
  next FDTD update. With identity rotation no coordinate transformation is needed.
- Native snapshot E(n-1) = -(Ac(n)-Ac(n-2))/(2dt), H derives from centered spatial
  differences of Ac(n-1). Native output snapshots are written before TDDFT updates
  the current; therefore snapshot J is the preceding-step charge current.

## Results

Compared all 257 FDTD snapshots (including initialization), every macro grid
point, all A/E/H components and preceding-step charge currents, plus 256 current
samples at each of the two TDDFT points.

| Comparison, max absolute error (a.u.) | PZ | TB-mBJ |
|---|---:|---:|
| A/c | 5.421010862427522e-19 | 5.421010862427522e-19 |
| E | 5.421010862427522e-19 | 5.421010862427522e-19 |
| H | 5.421010862427522e-19 | 5.421010862427522e-19 |
| Current (largest across both points) | 5.082197683525802e-20 | 5.082197683525802e-20 |
| Independent A-history replay current difference | 0 | 0 |

Replaying each recorded local A history into a fresh worker exactly reproduced
its coupled-run current, verifying independent local histories. A separate native
FDTD run with zero feedback current differed from the closed loop by up to
2.4990400307292674e-5 a.u. (PZ) and 2.5081863881438836e-5 a.u. (TB-mBJ) in A/c;
the feedback is nonzero, not merely a prescribed-field TDDFT test. Peak electronic
currents were 1.8529767904520924e-5 and 1.882890684854815e-5 a.u., respectively.
The final artifact script additionally checks the centered E identity and exact
one-step charge-current sign/lag.

## Build and limits

All application objects and the bridge passed the forbidden-eigensolver audit.
The Mac MPI C wrapper defaults to clang; `OMPI_CC=/opt/homebrew/bin/gcc-15` was
needed for the configured OpenMP flags. Local MPI execution requires permission
to open local communication sockets. Neither issue required physics-code changes.

SHA256:
- Native MPI reference: `9bfea204c38d3f4577695fdff52635274028c3e0e27e63ef955cca3af4d57ad0`
- Serial TDDFT worker: `c441bd0d922211b3fc48f47f594903e68811da539aba62cab0a36481a2ba8306`
- Native Weyl bridge: `f1331d173277633b6350a89eaa16c21b2b50b1e2a7436170e71c2ef51aafb5f8`

This is the existing A-form multiscale path. It is not validation of the standalone
E/H FDTD module coupled to TDDFT, nor a public E/H step interface. Native E/H
observables are nevertheless compared. Rotation, oblique incidence, 3D structures,
PML/EH coupling, time subcycling, restart of the entire coupled state, spatial/time
convergence and production scaling were not tested. The short/thin two-cell system
is a connection regression fixture, not a converged material simulation.

Inputs, outputs, binaries, bridge source, scripts, checksums and plots are saved at:
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/salmonpy-ms-two-points-20261113`.
The validation bridge is intentionally fixture-specific, not a supported public API.

Environment note: absolute paths in this report identify the historical validation host. They are not portable setup defaults or required user paths.
