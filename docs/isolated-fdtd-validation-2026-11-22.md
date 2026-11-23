# C2H2 isolated TDDFT / FDTD connection

User authorized the C2H2 campaign's GS, worker copies and outputs. Canonical
record: calc2026/salmonpy-c2h2-feedback-20261122. No new heavy CI run.

## Actual electronic reference

Existing audited executable SHA256
7f8e0c7f8d905413a8186969d3fe2191c6d583957d4c87d29eca26237c627017
was checked and reused unchanged. No source/object audit repeated at job startup.
PZ, fixed-ion C2H2,24×24×32,0.35Å, six states, ten electrons. One GS calculation:
123 iterations, residual1.0780440e-10<1e-9. GS saved for reuse.
Two20-step batch TDDFT references repeat exactly. Isolated E-worker/batch dipole
maximum difference3.85523175341039e-12 e·bohr. Repeated workers, interleaved
independent states and accepted-history replay match exactly.

## Two-field feedback

Four outer MPI ranks: one field rank plus bookkeeping rank per field context.
One fixed molecule, length-gauge endpoint E; full field minus own-source vacuum
reference. Synthetic LD is present only in the environment field. One initial
zero-current field interval; electronic time=field time-dt, reported explicitly.
64 field steps,dt=.003au. Fortran forms the charge-current source from dipoles;
Python coordinates trials. No change to existing TDDFT/FDTD physical update order.

Source vacuum and LD runs match independent electronic E-drive replay exactly.
Reject once at interval30, restore and retry: every retained E,p,j,drive and time
array equals continuous execution exactly. Installed-wheel vacuum/LD retry
histories equal source continuous histories exactly. Electronic restoration uses
accepted-history replay; direct/no-recomputation restoration is not claimed.
Field tokens are single-use; protect the accepted attempt with a fresh token.

The selected-E difference between LD and vacuum reaches1.14808e-5au. The own-source
vacuum E peak is1.00702e-6au; vacuum dipole change reaches7.92686e-6e·bohr.
These are nonzero wiring diagnostics, not converged response predictions.

## Actual-MPI defect found and fixed

The initial MPI-linked child worker inherited its parent's OMPI/PMI/PMIX/PRTE
connection information and waited during initialization. Initial test logs are
retained; only those test launchers were stopped. Other jobs were not modified.
The child's environment now omits those connection/rank entries and
MPI_LOCALRANKID/MPI_LOCALNRANKS/SLURM_PROCID/SLURM_LOCALID. Parent MPI/environment
and scheduler allocation IDs are unchanged. A synthetic contract test reproduces
the leak before the fix and passes afterward. Actual OpenMPI singleton workers
then complete in the four-rank fixture. Other MPI/site behavior is unvalidated.

## Scope

Source and installed-wheel suites:683 tests each, one optional skip each.
Paired guides and AI JSON are synchronized. No public general two-way driver,
quantum norm/work conservation, total electromagnetic energy certificate,
converged optical response, multinode or other-platform claim. No implicit SCF,
network download or build was introduced into runtime/import paths.
