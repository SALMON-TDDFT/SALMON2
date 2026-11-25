# Streamed E/H and electric spatial feedback validation — 2026-11-25

The public `IsolatedTDDFTState.step_grid_field` stages two spatial E/H frames for step n, using the existing integrator and direct checkpoint format. Geometry, dt, target index and complete frame buffers are checked before installation. The first nonzero magnetic input additionally checks the actual halo/projector query domain. A compact grid covering that domain is accepted. The tests do not claim a measured negative magnetic-coverage case.

`GridFieldSeries.difference` subtracts the own-vacuum reference in Fortran. The development tool `tools/test_grid_feedback_mpi.py` combines this with density-derived electric path currents and direct FDTD E-component injection. It uses two FDTD groups and one molecular worker, with the molecule displaced from FDTD grid points. Only field leaders transfer patch samples/current ownership; full grids are not broadcast across the application.

## Measured checks

- Final solver: Git-only export/preparation/build succeeded, including prohibited-eigensolver source/object audits. Source inventory SHA256: `9b747978eac5651e779c50b96695cb8d50fcff93e7156e19f838845ec3556d9d`.
- Executable SHA256: `44e3f6d8f8e18522aebed461e52de36aa9e0dc5cf38707cc9463e318fa6578c8`.
- Source-distribution-to-wheel build and clean installation passed. Wheel SHA256: `7216482ab37cb4ece1a24b61bbacbb0f24d2b62debbcb93239b57b9320e6ec66` (CPython 3.14, Apple Silicon/macOS 26).
- Installed targeted tests: 113 + 16 = 129 passed. These counts describe this turn's selected tests, not the entire test suite.
- Real C2H2: streamed/prepared E/H all nine observables are bit-exact. Zero-H negative control differs by 2.6330582159062033e-10 in the measured observable; malformed geometry/Infinity leave the state intact. Process-exit direct restart at molecular step2 to4 matches without extra evolution.
- Installed public `isolated_response --grid-series` sample agrees with the accepted-step reference. Its archive also includes the initial frame.
- Electric feedback: outer4/6-rank vacuum checks agree, including explicit rejection/retry and disk restart. The longer LD case uses192 molecular updates and193 FDTD updates, with a declared initial field hold. The returned environment field differs from illumination-only by about3.999e-7 au, above roundoff. Predicted/accepted drive disagreement is at most1.388e-17 au.
- Final-source and installed-wheel six-rank runs restore molecular96/FDTD97 with zero extra evolution. The installed run additionally rejects/retries step128. Final arrays, observables and all six ordinary `.data` streams are bit-exact with continuous execution.

## Scope and remaining work

This is electric-only feedback under an explicitly chosen own-vacuum subtraction policy. The one-way input accepts E/H, but independent magnetization-current return is not implemented. The development driver aborts unexpected MPI errors; it does not provide the general production coordinator's composed-failure guarantees. Discrete charge continuity, physical convergence, realistic tip/bow-tie geometry and rolling field-storage/I/O optimization remain separate work. The current full prepared field allocation scales with maximum steps and patch size.

The authorized existing C2H2 GS was reused; no new GS calculation or CI run was performed. Results reside in the canonical `calc2026/salmonpy-c2h2-nearfield-20261125` campaign: `grid-stream-wheel-live`, `grid-stream-wheel-example`, `grid-feedback-release-ld`, `grid-feedback-release-ld-restored`, and `grid-feedback-wheel-ld-restored`. Machine-readable record: `grid-field-feedback-validation.json`. Storage is local OneDrive; cloud synchronization was not independently confirmed. This is a software-connection/restart validation, not a physical-convergence claim.
