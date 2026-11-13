# TDDFT response-service connection, 2026-11-13

Implemented TDDFTMaterialProvider, TDDFTPreparedState and TDDFTResponseBackend in
the existing experimental distributed response service. The provider snapshots a
prepared GS lazily, detects changes since identity creation and retains immutable
initial data. Each macro point creates its own existing serial TDDFT worker.
The adapter supplies endpoint A/c and returns physical charge current (-native J).
No native TDDFT update routine or integration order was changed.

The backend independently verifies a managed runtime's tddft/single-rank
capabilities, matches its binary SHA to the initial-state descriptor and actually
launches that verified executable (not a different same-hash path supplied by a
custom loader). Custom providers can return TDDFTPreparedState with their snapshot
and managed runtime; raw unverified initial states are rejected.

## Local validation

Apple Silicon macOS, Python 3.14, Open MPI 5.0.9, mpi4py 4.1.2, existing audited
serial native runtime; one OpenMP/BLAS thread per worker. The GS fixture uses
synthetic normalized orbitals from test_tddft_entrypoint.py and is not a converged
material simulation.

- Three Python MPI ranks: one prescribed-field role and two one-rank material
  groups, two independent electronic states per material group (four total).
- Three different endpoint-A steps per point. All returned charge-current arrays
  are exactly equal to four separately evolved native reference states.
- A second session from the same provider reproduces every current exactly.
- Closing the initial-state provider after creating the second session leaves all
  live point states usable; they own their data/workers.
- Changed input bundle is rejected. A multi-rank material group is rejected before
  starting a serial-only backend. Unverified initial descriptors are rejected.
- A descriptor with matching binary SHA but a missing alternate executable path
  successfully uses the verified managed runtime; the pre-fix test failed.
- Final sdist -> repaired wheel -> isolated installation: 95 tests pass. From
  that installed wheel, custom-response MPI6, LD MPI4, synthetic example MPI5,
  TDDFT MPI3 and the shipped TDDFT example MPI3 all pass. MPI launchers are bounded
  and terminate their process groups on timeout.
- Review findings about serial-runtime validation and actual executable provenance
  were fixed, regression-tested and reviewed again with no remaining blocker for
  this stated scope.

CI configuration now invokes this integration against the existing synthetic
TDDFT relocation fixture on Linux. The new CI changes have not been pushed or run;
these results are local validation only.

## Scope

This is prescribed-A response exchange, not FDTD propagation. It uses a resident
serial subprocess for each macro point; MPI rank count is not total OS process
count. Plan memory for all states. Only one MPI rank per material group is
supported by this backend. Embedded multi-rank TDDFT and multiple electronic
contexts inside one Fortran process remain unimplemented. The next extraction
boundaries are recorded in plans/2026-11-13-tddft-context-inventory.md.

User sample: salmonpy/examples/mpi_tddft.py. Setup and replacement contracts:
salmonpy/DISTRIBUTED.md and DISTRIBUTED.en.md.

Historical local artifacts (not portable defaults):
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/tddft-response-service`.
