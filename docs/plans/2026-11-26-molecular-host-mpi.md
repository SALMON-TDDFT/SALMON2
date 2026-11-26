# Molecular two-level MPI (2026-11-26)

The caller starts one MPI application and splits its allocated ranks into a field
group and a molecular group. The molecular group supplies its communicator to
Fortran. No nested launcher or molecular child process is used.

## Implemented

- `run_molecular_mpi` runs isolated, fixed-ion, length-gauge molecular TDDFT on
  dedicated caller-owned ranks. `MolecularMPIController` supplies full spatial E/H
  frames. Commands and fields are distributed only within the molecular group;
  the caller retains MPI ownership after detach.
- The public `molecular_fdtd_mpi_response` example joins two FDTD ranks and two
  molecular orbital ranks. Coupling is one-way FDTD -> molecule, with no molecular
  current return or own-vacuum subtraction. Molecular step n uses field frames
  n,n+1. The molecular position may lie between field grid points.
- Direct checkpoints have rank-specific payloads and explicit decomposition,
  input and binary identities. Joint restart restores field/LD, molecular state,
  output prefixes and step counters without extra evolution.
- Shared packets are read once by the molecular root and broadcast within that
  group. A shared working filesystem is required. This is not an I/O-free channel.
- Bilingual guides, AI manual and an installed-wheel example document the route.

## Verification and limits

See `../molecular-host-mpi-validation-2026-11-26.md`. Control-only subgroup tests
were the first dependency; they were followed by actual C2H2 propagation and
process-exit joint restart. The earlier one-electronic-rank restriction still
applies to the separate-process worker, not this explicit host-MPI interface.

Validated scope: CPU PZ, fixed ions, spinless molecule, orbital decomposition
with one/two orbital ranks or spatial layouts (2,1,1), (1,2,1), (1,1,2),
and (2,1,2), plus two field ranks on one Mac. Dedicated main-thread
ranks support one molecular lifetime per process. Restart currently requires the
same decomposition and matching binary/input identity.

Spatial validation and same-layout direct restart are recorded in
`../molecular-spatial-mpi-validation-2026-11-26.md`.

Still open: combined molecular orbital/spatial decomposition, multiple nodes and larger
rank counts, several molecular states per process, density/magnetization exports
through this channel, and general failure-atomic recovery across field and molecule.
Unexpected distributed failures abort the example; no arbitrary MPI-failure retry
guarantee is implied. Physical convergence and bidirectional coupling remain separate.
