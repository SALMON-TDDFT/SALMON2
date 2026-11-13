# Two-level MPI response foundation validation — 2026-11-13

Scope: Python group orchestration, real SALMON communication-layer attachment,
independent user/LD macro-point states, prescribed-A distribution and physical-J
collection. This is not embedded TDDFT, an E/H FDTD coupling algorithm, layer
meshing or complete multi-material electronic dynamics.

Implemented interfaces are described in salmonpy/DISTRIBUTED.md and
DISTRIBUTED.en.md; the larger design remains in
plans/2026-11-13-two-level-mpi-interface.md. Full-solver modifications remain in
the pinned external runtime patch. No upstream source tree was restored.

## Verified locally

macOS arm64, Python 3.14, mpi4py 4.1.2, Open MPI 5.0.9, GNU Fortran 15.
One OpenMP/BLAS thread per rank. Memory was not measured; these are tiny regression
fixtures, not throughput/scaling benchmarks.

- Native MPI, five ranks: two-/three-rank groups have independent rank/size,
  reductions and unequal numbers of default barriers. Attach duplicates the
  subgroup; attached legacy comm_init does not initialize MPI twice;
  comm_finalize detaches without finalizing the MPI owner. Duplicate attachment
  and detached-state errors are rejected.
- Python/Fortran, six ranks: one field-role rank plus two-/three-rank material
  groups. Four independent user states share an immutable initial definition.
  A/J histories match the manufactured serial law. Regressions cover inconsistent
  plans, duplicate material-key definitions, strided result buffers, collective
  reduction overflow, failed-session rejection, successful-state cleanup after
  factory failure, duplicate native attachment and premature context closure.
- Shared-lifetime regression: native contexts exist only on material groups;
  every group rejects context closure while any group is busy. Pre-fix test failed;
  after the fix no leader communicator is freed prematurely. Busy status is agreed
  across leaders and broadcast within every group before any free.
- Distributed LD, four ranks: two field ranks plus two one-rank LD groups;
  different pole counts, two independent points per material. Five endpoint-A
  steps match sequential LD current arrays exactly. This adapter uses backward
  difference E and is not the legacy implicit FDTD LD/Q update.
- Shipped user-code example, five ranks: all four response histories written to
  a caller-selected JSON file. Optional Fortran attachment also executes.
- Full MPI SALMON rebuilt from the updated external patch: all 19 existing FDTD
  material-interface comparisons pass, including old/new singleton pulse/impulse
  outputs and orientation/ID permutations. Standalone numerical order is retained.
- Checkout Python suite: 94 tests pass. Fresh sdist -> repaired wheel -> isolated
  installation: 94 tests pass. Installed-wheel tools successfully prepare the
  external source and build an audited host-MPI bridge. The installed wheel with
  that bridge passes the six-rank, four-rank and five-rank tests again.
- Shared LD/SBE kernel equality and prepared patch checksums checked. Full native
  application and bridge source/object forbidden-eigensolver audits pass.
- Review blockers fixed and re-reviewed: material cache identity, contiguous MPI
  buffers/nonfinite reductions, constructor cleanup, bridge audits and shared
  communicator lifetime. Collective backend factories must succeed on all group
  ranks or clean partial resources and fail collectively; this is an explicit
  user-backend contract.

Runtime-bundle CI is configured to build/test the foundation on Linux with
OpenMPI, including time-bounded launchers and explicit LD kernel preparation.
These new CI changes have not been pushed/executed. Local success does not claim
cross-platform MPI validation. Ordinary LD wheel compatibility remains separate.

## Reproduction

Build the explicit bridge from newly prepared runtime source, then run:

```sh
PYTHONPATH=. python tools/run_parallel_mpi_tests.py BRIDGE_BUILD_JSON NEW_RESULTS
```

The runner uses six, four and five ranks and terminates its process group on
timeout. MPI must be installed consistently for Python and Fortran. It does not
compile on import or download an MPI runtime. See the paired guide for setup.

Historical local artifacts, not environment requirements:
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/two-level-mpi`.
Artifacts include sample JSON, source/build provenance and local verification
logs. Checkpoints, physical Maxwell-current coupling and production material
validation are not provided by the synthetic demonstration.
