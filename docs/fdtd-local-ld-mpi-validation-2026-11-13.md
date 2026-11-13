# Existing MPI FDTD with rank-local LD, 2026-11-13

User-approved allocation: electronic solver groups initially own one material and
one macro point; LD-dominated calculations instead partition the large FDTD grid
and keep many independent LD histories on each field rank. No dedicated LD MPI
group is allocated per cell. SALMON's existing MPI Maxwell and implicit E/LD/H
numerical updates are used unchanged by this launcher.

Implemented FDTDMPI(grid, launcher), SalmonExecutable.from_mpi_build(build.json),
and the mpi= argument to the FDTD batch runner. The completed build record stores
MPI capability alongside binary SHA; later CMakeCache changes cannot turn a
serial build into an accepted MPI build. Portable runtime packaging similarly
requires a recorded serial build. Trusted external binaries may use the existing
low-level constructor with an explicit mpi_enabled=True capability assertion.

The controller must be a fresh serial Python process. Nested rank environments
and already initialized mpi4py are rejected, including direct execution of an
MPI-enabled binary without mpiexec. Interruption/timeout cleans the launcher
process group and records cancellation/failure rather than leaving started status.
Remote scheduler cleanup remains the launcher's responsibility.

## Local verification

macOS arm64, GNU Fortran/OpenMPI; one OpenMP/BLAS thread per rank. New MPI and
serial full-runtime builds completed source/object audits and correctly recorded
mpi_enabled true and false. The serial runtime also packaged successfully and
ran the revised TDDFT example.

- Three synthetic LD layers with 1, 2 and 3 poles; grid 48x40x40, spacing 2 a.u.,
  dt=0.003 a.u., 160 steps. Each material occupies 256 cells (768 matter cells
  total, 76032 vacuum cells). This is finite synthetic geometry, not a converged
  film or production-scaling benchmark; memory was not measured.
- Existing native FDTD on one rank versus four ranks with process grids (2,2,1)
  and (1,2,2): all three observer E/H arrays are exactly equal.
- Reversing material IDs and material tuple consistently on four ranks preserves
  every observed array exactly.
- A real MPI run with a short timeout raises TimeoutExpired and records failed;
  the launcher exits. Unit regression confirms interrupted runs record cancelled.
- Build-capability regression rejects a serial completed record even if the
  current cache is changed to MPI ON. Changed executable SHA is rejected.
- 100 Python tests pass from checkout and from a fresh sdist -> repaired wheel ->
  isolated install. Installed-wheel MPI FDTD regression and three-material sample
  both pass. The wheel's standard TDDFT example runs with two material groups,
  one material/one point per group (two total electronic states).
- Review findings on build provenance, nested MPI and cancellation were fixed
  and re-reviewed without remaining blockers for this launcher scope.

CI is configured to run the MPI FDTD regression on Linux. These current changes
have not been pushed/executed on CI; local success is not cross-platform proof.
MPI launching currently targets macOS/Linux. LD wheel platform support is separate.

Usage: salmonpy/FDTD_MATERIALS.md and FDTD_MATERIALS.en.md.
Sample: salmonpy/examples/fdtd_mpi_layers.py (start from serial Python).
Electronic sample: salmonpy/examples/mpi_tddft.py (start with three Python ranks).
The generic response service retains its tested multiple-state capability, but
current electronic examples use the user's simpler one-point-per-group profile.
Full E/H FDTD/TDDFT coupling and inner-MPI electronic contexts remain pending.

Historical artifacts, not environment requirements:
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/fdtd-local-ld-mpi`.

## Documentation and sample follow-up

Japanese/English onboarding and sample guides now distinguish the serial-Python
MPI-FDTD launcher, collective user-response example, and serial-worker TDDFT
response service. They state rank/state ownership, required runtime types,
postprocessing outputs, resource limits and unavailable full coupling. AI task
`distributed_tddft` separates its runtime requirements from user callbacks;
`sample_files` and embedded guide snapshots are synchronized.

The custom response sample now uses one point per two-rank material group
(point IDs 10 and 20) and two different material coefficients. The MPI regression
first failed against the former four-point sample, then passed with the new
sample. Every time frame is compared with its analytic current reference.

A fresh source distribution, repaired macOS arm64 wheel and isolated installation
passed 100 Python tests. Installed-package custom/LD/TDDFT MPI regressions and the
five-rank custom, three-rank TDDFT and four-rank FDTD examples passed. The FDTD
example retained 256 cells for each of its three materials. Wheel documentation
and sample bytes match the checkout. This follow-up changes guidance and the
manufactured sample; it does not implement embedded TDDFT or E/H coupling.
