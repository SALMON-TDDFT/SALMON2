# Two-level MPI response foundation (experimental)

[日本語](DISTRIBUTED.md)

Python splits the launched MPI job into solver groups. A material group owns one
or more stable macro-point IDs. Each point gets an independent response state.
The field group distributes prescribed A/c and receives physical charge-current
density in point-ID order. All values use atomic units; A and current correspond
to successive fixed-dt endpoints, with initial A=0. This is a response service,
**a prescribed-A service by itself**. Add the [solver coupling entries](COUPLED.en.md) for internally parallel TDDFT and E/H field stepping.

Install the optional `salmonpy[mpi]` dependency in an environment with a compatible
MPI installation. Ordinary import of SALMONpy/distributed does not import MPI,
compile code or start ranks. Start the Python job once with mpiexec or your
scheduler; each MPI rank runs Python and calls the same collective interface.

```sh
mpiexec -n 5 python -m salmonpy.examples.mpi_response NEW_JSON
```

The sample uses one field-role rank and two material groups of two ranks each,
with one macro point per material group (IDs 10 and 20) and different material coefficients. Its user-defined memory
law is synthetic; it is not a physical material model or FDTD calculation.

## User-supplied material preparation and response

`MaterialBinding(key, provider, backend)` separates the two replaceable objects.
Reuse the same binding instance for points using one material key. Provider
`identity` must identify its version **and configuration**, consistently on every
rank. `prepare(context)` returns provider-owned immutable initial data. It may
load any user format; no SALMON data structure is required. calculation implementation `identity`
has the same consistency requirement. calculation implementation implements:

- `execution`: `collective` (all group ranks execute) or `leader_only` (currently
  requires a one-rank material group).
- `reduction`: `sum` (one local contribution per rank) or `leader` (leader already
  holds the complete current; other ranks still return finite three-vectors).
- `create_state(initial, point_id, dt, context)`: returns a fresh state with
  `advance(A_old, A_new)` and `close()`. advance returns a finite three-vector of
  **physical charge-current density**, not electron-number current or total current.

Callbacks receive `SolverContext`: only their solver communicator, group name and
point IDs. All members call the same points in the same order. Do not access
MPI_COMM_WORLD for calculation implementation collectives, initialize/finalize MPI, mutate shared
initial data, or return one mutable state for several points. State close is
collective within its group when the implementation requires it. A failing
prepare/create callback must collectively release its own partially constructed
resources before raising. The service releases previously completed states.

`ParallelPlan` validates exclusive rank/point ownership and exactly one field
group. `ParallelContext.create(world, plan)` collectively checks plan agreement
before Split. `ResponseSession(context, bindings, dt=...)` collectively prepares
states. Every rank calls `advance`; only the field-group leader supplies a
float64-compatible `(number_of_points,3)` array, others pass None. Every field
rank receives J in `plan.point_ids` order; material ranks receive None.
Providers/calculation implementations must not change their identity/configuration while a session
is active. Contexts, communicators and states are not thread-shareable.

Invalid drive inputs are rejected before evolution. A calculation implementation step failure marks
the session unusable on every group; there is no cross-group rollback. Close the
session and restart from a known initial state. solver process failure or a user
callback deadlocking inside MPI requires launcher-level timeout/cancellation;
Python exception agreement cannot recover it. Close sessions and detach from calculation libraries
before closing ParallelContext. MPI finalization belongs to the Python job owner.

`LDMaterialProvider(model)` and `LDResponseBackend()` reuse the existing LD state.
They compute E_new=-(A_new-A_old)/dt and return physical J from LDState.step.
Their identity is `salmonpy.ld.a-endpoints/1`; this **does not implement the
semi-implicit legacy FDTD LD coefficients/Q coupling**. Do not substitute it into
that update schedule. The public TDDFTState/SBEState current is already charge
current j; do not negate it again. A raw J_matter observable needs one explicitly
declared conversion. The response service never guesses the meaning of user current.

## Fortran subgroup attachment

Build explicitly from a newly prepared external source, using the same MPI
installation as mpi4py (currently GNU-compatible mpifort, macOS/Linux):

```sh
python tools/prepare_salmon_runtime.py NEW_SOURCE
python tools/build_parallel_bridge.py NEW_SOURCE NEW_BRIDGE --fc mpifort
```

An installed wheel includes these tools under `salmonpy/_runtime_source/tools`;
locate the directory with Python's importlib.resources. Builds keep source/object
audits and write build.json, hashes and logs. The MPI bridge is calling program-specific and
is not bundled with a serial runtime or automatically compiled. Python and
Fortran must link a compatible MPI ABI; a runtime-name/initialized-state check is
an additional guard, not universal ABI certification.

Pass its library path and build.json SHA256 to
`NativeParallelContext(context, library, sha256)`, or use the sample's optional
`--bridge PATH --sha256 HASH`. The bridge accepts the process-local Fortran handle
from mpi4py, checks integer representation, and the SALMON communication module
duplicates the supplied subgroup. Default communication and barrier then use
that subgroup. Detach frees only the duplicate, never MPI itself. `info`, `sum`
and `barrier` expose attachment diagnostics; they are not solver entry points.
One solver attachment per process/library is currently allowed. The solver
MPI test also checks legacy comm_init/comm_finalize while attached.

One-point TDDFT called in the Python process and restricted E/H↔A/J coupling are available through
[COUPLED.en.md](COUPLED.en.md). Remaining work includes in-process DFT, safe multiple
electronic states per process, automatic layer meshing, restart and production scaling. Existing batch FDTD and process-
isolated TDDFT APIs continue to serve their documented scopes.

For a collective create_state callback, every group rank must either return successfully, or all must release partial resources and raise. Mixed success/failure across ranks is outside the callback contract. Context close is collective over all groups and checks that no group still has dependents before freeing the shared leader communicator.


## Process-isolated TDDFT adapter

`TDDFTMaterialProvider(bundle, input_file, runtime=managed_runtime, max_steps=...,
threads=1)` verifies the managed runtime's `tddft` and `single-rank` capabilities.
Its identity includes input-file contents, assets, binary checksum and worker
settings. Preparation lazily snapshots the bundle in each material group and
rejects changes since identity creation. Close the provider after session use;
already-created states own copies and survive provider closure. Each point gets
its own existing TDDFT worker, including wavefunctions, Hartree history and output
directory. No new propagation algorithm is introduced.

`TDDFTResponseBackend()` is `leader_only`: **each material group must have exactly
one MPI rank**. One group may own several points, processed sequentially. Every
point retains a separate resident serial subprocess and its memory/files. Budget
memory for every point and CPU threads for active workers; the MPI rank count is
not the count of all OS processes. This serial-worker calculation implementation does not provide
internal MPI; use the [separate solver calculation implementation](COUPLED.en.md) for that and E/H coupling. `TDDFTState` already
returns the physical charge-current density j = −J_native; the service does not negate again.
The supported XC/propagator scope is exactly that of the existing TDDFT worker.
It loads a prepared GS; it does not run DFT preparation automatically.

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_tddft RUNTIME GS_BUNDLE RT_INPUT NEW_JSON
```

The standard example keeps two independent electronic states across two one-rank material
groups. It supplies short synthetic A ramps, not a validated laser experiment.
Converge the GS, spatial grid, k points, bands and timestep for real calculations.
The example aborts its entire MPI job on an error to avoid stranded peer ranks.

Custom loaders returning an existing TDDFTInitialState must wrap it in `TDDFTPreparedState(initial, managed_runtime)`. The TDDFT calculation implementation independently rechecks serial-runtime capability and executable checksum before spawning. An unwrapped or mismatched initial state is rejected.

Current working configuration: one material and one macro point per electronic MPI group. The API retains independent multi-state capability for future work, but the standard TDDFT sample assigns one point per group. LD-dominated Maxwell calculations instead use solver spatial MPI FDTD with many local LD states; see FDTD_MATERIALS.en.md. These are distinct execution paths.

For one serial model contract, the explicit conformance checker and the point-MPI adapter, see [MATERIALS.en.md](MATERIALS.en.md).

For explicit identity propagation/output in MPI response sessions see
[MPI_IDENTITY.en.md](MPI_IDENTITY.en.md). PointRun itself remains serial; restart is not provided.

Before provider preparation, ResponseSession validates optional calculation implementation.supported_group_sizes: a nonempty tuple of unique positive integers. Declarations must agree across all groups; tuple order is irrelevant. Undeclared calculation implementations retain existing execution/reduction checks, without capability certification. The A-packet wrapper forwards the original declaration.


[Shared HSE multipoint placement and protected direct resume](SHARED_HSE.en.md#normal-mixed-mpi-sessions).

For opt-in standard TDDFT shared models and multipoint groups, see [SHARED_TDDFT.en.md](SHARED_TDDFT.en.md).


## Minimize outer collectives

Application-wide communication is limited to necessary setup, consistent failure
agreement and explicitly requested final aggregation. Runtime validation record reuse is
rank-local. The validation-only mpi_runtime_receipt example uses one setup bcast
and one final gather to rank0, without allgather, health allreduce or Barrier.
The final gather itself shows caller MPI remains usable; reports are not copied
to every rank. The script terminates after reporting and never starts a solver.

NativeEHField.advance coalesces its Python wrapper's subgroup checks from six
collectives to two: (local input error, current values), then (solver status,
finite E, finite candidate A, integration error). Every member sees the same
records and decides identically. No reduction of solver physical communication
is claimed. Signed-zero equality, before-entry rejection, recoverable solver
rollback and commit-after-agreement are preserved. No numerical update changes.

Keep model-internal physical sums/halos and required failure decisions.
ResponseSession now scatters A rows only to their owning material leader, then
broadcasts those rows only within that material group. Field groups receive no
response-service input payload; singleton groups need no subgroup broadcast.
Field packets are similarly scattered as per-owner dictionaries, preserving
values, units, positions and clocks. J already gathers local owner contributions;
the callback and solver update order is unchanged. solver spatial return routing is described in COUPLED.en.md.

Root need not be first in plan order; point IDs and owning groups may interleave
and own different counts. Cached routing indexes restore global J order without
repeated per-step point searches. communication_report records received_A_rows_by_group;
logical_A_bytes_per_step is the logical dataset size, not measured network traffic.
No new world collective is added. Existing error/finite agreements are retained.

For checkpoint compatibility response._A retains a global-indexed shape: only
owned rows are current on material ranks, and only the field leader retains full
response-service A history. Nonowned/private rows are placeholders, not replicated
output. Use field-leader point archives for global observations. Old checkpoints
with full rows remain readable. solver field A/E and physical state data are not
replaced by placeholders. solver spatial J routing is described in COUPLED.en.md.
Adjacent outer input checking is also coalesced as described in COUPLED.en.md; full solver communication is not claimed optimal.


## Spatially owned FDTD called directly from Python current routing

With a matching runtime, CoupledStepper and PacketCoupledStepper query
`sp_eh_current_mask` once after field initialization. A point current is sent by
Scatterv to every field rank owning an edge touched by its transpose sample-map
deposit, rather than only the rank owning the point center. Boundary-crossing
points reach both sides; empty owners receive zero current entries. The three
components of an owned point travel together; this is not component-level compression.

The field leader retains globally ordered J for output and checkpoints. Other
field ranks retain the global array shape, but only owned rows are meaningful.
Use leader output for full-point analysis. A/E samples, edge weights, the
E→LD→H update order, collective failure decisions and restore guarantees are unchanged.

The new `sp_eh_step_owned` ABI accepts packed local point currents. The existing
`NativeEHField.advance(J)` still requires identical full-point J on every field
rank. Older runtimes without this ABI, Multiscale1D and standalone ResponseSession
keep replicated results. Build explicitly; importing or starting a calculation
does not rebuild a runtime automatically.

Checkpoint formats and MPI compatibility rules are unchanged. Full-row J saved
with the same compatible runtime and decomposition can resume by selecting the
owned rows. This does not bypass binary or MPI-layout compatibility checks.
The field manifest records `current_routing` and `current_rows_by_field_rank`;
explicit `communication_report()` rank records expose `received_current_rows`.
These are logical point counts, not measured network bytes.

See COUPLED.en.md for status-bearing A/J transport, failure-only explanatory Bcast and current communication counts.
