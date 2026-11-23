# Single-state TDDFT performance mode — 2026-11-23

The explicit NativeTDDFTBackend execution_mode selects the single-state host-MPI
route. Its Fortran time-evolution call is shared with checked mode; the checked
path retains rollback/diagnostics. The performance path omits backups, finite
scans and diagnostic agreement, retaining setup checks, physical solver work and
explicit observation/restore. Returned errors mark the Python state failed;
unknown Fortran/MPI failures may terminate the job. Default remains checked.
Shared/registry models, fast HSE, process workers and molecules remain pending.
Outer FDTD/coupling validation is independent. No new numerical Python algorithm.

A small existing manufactured GS (8³ grid, 2³ k grid, two orbitals) was read via
symlinks. No GS calculation, GS copies or electronic-state disk checkpoints.
Three independent MPI groups compared old-library checked, new-library checked
and new-library performance over 16 steps. Explicit in-memory capture at step8,
restore and continuation gave identical currents without replaying the prefix.
PZ/TB-mBJ/Libxc-PZ/PBE × 1/2 ranks per group passed from source and installed
wheel (16 cases). Current and measured electron number matched exactly across
modes/binaries; installed/source arrays matched too. This is update regression,
not a converged material study, norm-accuracy recommendation or speed benchmark.

Git-only source preparation/new runtime build and source/application-object
prohibition audits passed. Source and installed-wheel suites: 701 tests each,
700 passed, one optional skip. Wheel code, paired guides, AI JSON and sample
match source bytes. Required CI was not dispatched manually.

Reproduce with tools/test_tddft_performance_mpi.py using an existing normalized
small fixture; see its module instructions. Records and artifacts reside at
calc2026/salmonpy-tddft-performance-20261123. Local OneDrive storage does not
confirm cloud synchronization.

Runtime SHA256: `d95d69e9ac06d47a6f15e92c5f228e1bd8c8eb228aa70b3104bd1d0b98bb485b`.
Wheel SHA256: `1c82ff91ef7c14db9b6249900b0ef8031eece0072542bc493c4a30fa7f630943`.
