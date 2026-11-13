# Two-level MPI material-response interface

Status: communicator/response foundation implemented experimentally; full coupled
solver design remains pending. See salmonpy/DISTRIBUTED.md and DISTRIBUTED.en.md
for actual available APIs. The proposed signatures below are design targets. Incorporates the user's explicit
requirements: layer order/thickness and resource budget; Python orchestration;
Fortran groups with internal MPI; one or more independent macro points per group;
A distribution and spatial J collection. No upstream merge is planned.

## Communicator boundary (first implementation milestone)

An MPI launcher starts the complete Python job once. MPI initialization creates
MPI_COMM_WORLD; Python does not manufacture MPI handles. A serial convenience
launcher may invoke mpiexec/scheduler exactly once. An already launched MPI job
must not recursively launch another job. Dynamic MPI_Comm_spawn is not required.

The proposed `salmonpy.distributed` module has two phases:

1. `ParallelPlan.validate(world_size, backend_capabilities)` checks the entire
   assignment before any split, on every rank. Return the same structured failure
   on every rank if assignments, material definitions or plan hashes disagree.
2. `ParallelContext.create(world, plan)` is collective over world. It calls
   `world.Split(color=group_id, key=world_rank)` for solver groups, and creates a
   separate leader communicator for FDTD/material exchange (nonleaders use
   MPI_UNDEFINED). All ranks execute splits in the same order.

Each group receives `solver_comm`, its local rank/size, role, stable group ID,
ordered macro-point IDs and a separate orchestration communicator. A Fortran
backend collectively duplicates solver_comm during attachment. Its duplicate is
`comm_solver`, the user's conceptual MPI_COMM_WORLD2. This is an ordinary MPI
communicator, not a new predefined MPI world or a second MPI initialization.
Fortran's existing k/orbital/real-space decomposition subdivides comm_solver.
Python owns the outer distribution, the backend owns its internal subdivisions.

Example resource plan: 12 total ranks = FDTD group of 4 + two TDDFT groups of 4.
Each TDDFT group can own several macro points; this does not multiply the number
of ranks again. LD may use one-rank groups or share a material group's sequential
point schedule. Report per-group ranks, threads/rank, number of states and memory
estimate/measurement/unknown status before execution. Reject unaccounted ranks or
an inner decomposition whose product differs from its group's size. Initial
support uses fixed assignments and no process migration or rank sharing between
simultaneously active groups.

### Python–Fortran ABI and lifetime

Proposed entry points (not implemented):

| Entry | Collective scope | Contract |
|---|---|---|
| `backend.attach(solver_comm, options)` | all ranks of solver_comm | MPI must already be initialized; duplicate communicator, query local rank/size, create backend context |
| `backend.prepare(context, material_spec)` | solver group | DFT or loader prepares immutable initial state |
| `backend.new_state(context, initial, point_id, step_config)` | solver group | create independent distributed mutable state; return opaque local handle |
| `backend.advance(context, state, drive_frame)` | solver group | same point and step on all members; return group-reduced response on leader |
| `backend.close_state(context, state)` | solver group | release state resources in identical point order |
| `backend.detach(context)` | solver group | wait for requests, free owned subcommunicators and duplicate only |

Use an MPI-aware compiled bridge, not an untyped ctypes integer cast. With current
SALMON `use mpi` integer bindings, mpi4py `Comm.py2f()` gives a Fortran communicator
handle. A compiled boundary must check/convert the MPI_Fint and Fortran INTEGER
representation; do not assume Python int, C int, MPI_Comm and mpi_f08 types share
an ABI. A C bridge can accept MPI_Fint and use MPI_Comm_f2c when a C handle is
needed. Any future mpi_f08 bridge needs its own supported conversion. Handles
are process-local and must never be serialized in plans/checkpoints or sent to
another process. Python and native modules must use a compatible MPI runtime
and ABI in the same process; reject mismatches before attachment. Serial runtime
bundles are not usable for this embedded MPI path.

Python/launcher owns MPI initialization/finalization. Fortran attachment does not
call MPI_Init or MPI_Finalize. Query initialized/finalized state and thread level;
initial mode confines all MPI calls to the main thread (no requirement for
MPI_THREAD_MULTIPLE). Dup/free are collective. Detach native contexts before
freeing Python communicators and before finalization. Imported Python modules do
not create communicators or compile anything.

Current source audit: external `src/parallel/parallelization.f90` setup_parallel
calls comm_init and comm_get_globalinfo; end_parallel calls comm_finalize.
`src/parallel/communication.f90` comm_get_globalinfo returns MPI_COMM_WORLD and
comm_sync_all without an explicit argument barriers on MPI_COMM_WORLD. Merely
passing a communicator to the worker is therefore insufficient. Embedded entry
points must attach an explicit context; all default communication, timing,
root-only output, reductions, BLACS/ScaLAPACK contexts and error paths must resolve
to that context. Retain existing standalone entry points and initialization order
for CLI compatibility. Audit all MPI_COMM_WORLD/SELF, MPI_Init/Finalize, STOP and
MPI_Abort uses, including wrappers, before enabling embedded capability.

## Material and macro-point state boundary

Separate immutable `MaterialSpec` (backend, coefficients or DFT input, XC/potential
and provenance), immutable `InitialState` (ground-state/loader result), and mutable
`MacroState(point_id, initial_id, state_handle, step, A_previous, histories)`.
A MaterialSpec ID is only a lookup key. Never use its numerical magnitude for
physics, ownership priority or rank assignment.

`LayerStack` specifies ordered layers and thicknesses with units. A discretization
step produces a MacroMesh: stable point IDs, physical coordinates, material IDs,
cell volume/weight and the mapping to the Maxwell grid. Point IDs are distinct
from material IDs, grid indices and MPI ranks. Report how thickness is rounded to
the selected grid; do not silently alter it. Vacuum points can omit a solver and
supply zero material current. Surface currents and vacuum propagation remain the
Maxwell solver's responsibility.

DFT prepares electronic initial states; TDDFT evolves them. SBE is another
A-driven backend. LD remains internally E-driven. Two points referencing the same
initial state must have separate wavefunctions/density matrices/P/J/Q, integration
history, observables and output directories. Immutable data may be shared only
where the native implementation guarantees it cannot be mutated.

A group owning several points executes them in deterministic point-ID order;
all ranks participate in each point's collective native call. One rank must not
independently skip to another point. Supporting multiple TDDFT states in one
process requires moving mutable SALMON module data into explicit contexts (or a
complete, tested save/restore boundary). A communicator split alone does not make
SALMON reentrant. Current TDDFTState uses one subprocess per state and is not an
implementation of this embedded collective interface. Advertise
`multiple_states_per_group=false` until isolation tests pass; reject incompatible
plans rather than share state silently. Achieving true multiple-state support is
a required milestone for the user's complete design, not an optional promise.

## A/J exchange contract

Use fixed binary float64 buffers, point IDs and explicit metadata. Proposed
`DriveFrame` fields: protocol version, run ID, coupling scheme, step, phase,
ordered point IDs, t_old/t_new, A_old/A_new and optional scheme-defined stages.
A is SALMON Ac=A/c in atomic units. Shape is (number of points, 3), with declared
layout. The old endpoint must match the state's accepted A. Reject duplicate,
missing, stale, out-of-order, nonfinite or unknown point/step messages.

`ResponseFrame` carries the same run/point/step/phase identity, physical charge
current density J_charge, its actual time position, status and selected
observables. The adapter converts native electronic J to charge current exactly
once (current SBE/TDDFT API uses J_charge=-J). Never return total current where
Maxwell expects density; cell-volume weighting belongs to the mesh transfer.
Do not assume LD and electronic backends share the same raw sign convention.

FDTD group computes/drives its distributed fields and uses a fixed transfer map
to sample A at macro points. Its leader initially sends batches to material-group
leaders; each leader broadcasts its batch inside the solver group. The backend
performs its native reductions internally and returns one response per point,
not one duplicate contribution per rank. Material leaders return batches, FDTD
validates completeness and inserts spatial current using a documented reverse
transfer map. A leader-mediated exchange is the first correct implementation;
later direct distributed routing must preserve exactly the same point contract.

Maxwell and macro grids need not coincide. The initial supported map should be
explicitly one-to-one or the existing native multiscale map. Interpolation,
volume weights, component locations and reverse deposition must be recorded;
a new averaging/interpolation rule cannot be inferred from layer thickness.
The E/H-to-A adapter owns A history, gauge/boundary convention, staggering and
sampling/deposition. Crystal backends receive only A, never an alternative E
forcing option. At this unit convention E=-d(A/c)/dt, but a generic finite
difference must not replace the existing coupling scheme's discrete rule.

## Preserve the numerical schedule

A single `advance(A_new)->J` is sufficient as the high-level material interface,
but is not enough to reproduce the existing semi-implicit LD/FDTD update.
LD contributes coefficients to E before E is updated and also maintains Q.
Do not replace that coupling with explicit J injection.

Declare capabilities and a fixed `coupling_scheme` at construction:

- Electronic A-driven scheme: the native multiscale scheduler determines which
  A endpoints/stages are sent and the time position at which J is consumed.
- `legacy_eh_ld`: prepare LD coefficients/old-Q contribution before E, then follow
  the existing old-J/E save → E update → source → old-Q injection → E exchange →
  LD J/P/Q update → linear response → old-H save → H update → H exchange sequence.
  Expose an additional internal `prepare_field_update` payload for coefficients
  and auxiliary sources. Its named fields, staggering and units must be mapped
  directly to the native LD kernel, not confused with physical J_charge.

Do not assign a numerical time position or convert the two schemes until native
fixtures establish the exact equivalence. A drive frame may contain additional
scheme-defined stages; subcycling/predictor-corrector iterations require explicit
iteration IDs and state semantics, not repeated mutating calls to `step`.
First coupled release uses fixed dt and a single validated scheme per run. Mixed
LD/TDDFT in E/H FDTD must pass its own native-equivalence tests before capability
is advertised.

## Failure, outputs and restart

All ranks of a solver group agree on precondition failures before entering native
collectives. Leaders report status to the orchestrator before the next global
step. A failed step must not be reported as accepted by another group. Initially,
a native failure stops the run and resumes only from a complete accepted-step
checkpoint; do not promise rollback without snapshots of every participating
state and Maxwell field. MPI process failure is not recoverable by ordinary
Python exceptions; timeout/job cancellation must terminate the entire job and
record the failure rather than leave other groups waiting.

Output paths are keyed by run/group/point, never rank alone. Store macro mesh,
material provenance, plan, A/J times and units, transfer map, communicator sizes,
and backend versions. Large densities, microscopic currents and wavefunctions
remain opt-in. Restart records all points and Maxwell state at one accepted step;
rank remapping requires a separately validated loader. No MPI handles in files.

## Acceptance gates and implementation order

1. Communicator attachment bridge: split two unequal groups; local rank/size,
   reductions and barrier tests; detach one group while the other continues;
   MPI still initialized afterward. Timeout detects accidental world collectives.
2. Stateful LD backend: two groups, several points each, independent histories;
   compare to sequential reference; permute point/material IDs and group mapping;
   A-to-E and physical-current sign/units tests.
3. TDDFT attachment: one point per group, then two points per group; same initial
   state with different A and repeated runs. Compare with current process-isolated
   worker and native inner MPI decompositions. Verify XC, potentials and all
   mutable state remain independent; retain source/object eigensolver audit.
4. FDTD collective context and phased coupling: retain singleton byte equality
   where operation order is unchanged, use stated tolerances for changed parallel
   reductions, then two-material native reference. Test A/E/J timing, transfer
   weights and LD implicit coefficients, not just final waveforms.
5. Layer-driven public workflow: dry-run resource/mesh report, paired human docs,
   AI JSON examples explicitly marked available only after implementation, and
   small end-to-end MPI CI plus opt-in full TDDFT integration.

Official API references:
- https://mpi4py.readthedocs.io/en/stable/reference/mpi4py.MPI.Comm.html
  (`Split`, `Dup`, `py2f`, `Free`).
- https://docs.open-mpi.org/en/main/man-openmpi/man3/MPI_Comm_split.3.html
  (collective subgroup construction and MPI_UNDEFINED).

## Module responsibilities for the first patch

Proposed files/entry points, deliberately not exported as usable APIs yet:

- `salmonpy/distributed.py`: `ParallelPlan`, `ParallelContext.create/close`;
  group/leader split, resource validation and collective lifetime.
- Native `salmonpy_parallel_context` module in the external compatibility patch:
  `attach_parallel(comm_f, context, ierr)` and `detach_parallel(context, ierr)`;
  MPI_Initialized/Finalized checks, communicator duplication, local rank/size,
  context-bound timing/error communication. Never reinitialize/finalize MPI.
- Native backend context: calls existing initialization through the supplied
  parallel context. All derived decomposition communicators originate there.
- `salmonpy/multiscale.py`: mesh/point ownership and A/J exchange orchestration;
  no electronic physics or implicit LD coefficient arithmetic here.
- Backend adapters: material preparation, independent point states, numerical
  schedule capabilities, physical-current normalization and observable selection.

Conceptual call sequence (design only):

```text
launcher -> MPI Python processes
  all ranks: validate plan -> split world -> create leader communicator
  each solver group, all members:
    Fortran attach_parallel(solver_comm.py2f()) -> duplicate comm_solver
    backend prepare -> new_state for each assigned point
  every accepted time step:
    FDTD group: prepare fields and drive payload at scheme-defined stage
    leaders: distribute A batches; each group broadcasts within solver_comm
    material groups: collectively advance each point in identical local order
    leaders: return point-keyed current batches
    FDTD group: validate and deposit J at scheme-defined stage, finish step
  each solver group: close states -> detach_parallel
  all ranks: free leader/group contexts -> MPI owner finalizes
```

`comm_solver` is the module-owned MPI_COMM_WORLD2-equivalent requested by the
user. The module runs on **every member** of the group, not only its leader.
MPI_Comm_dup does not create processes or increase parallelism. A backend can
have more internal communicators (k/orbital/grid), but cannot reach outside its
parent group unless explicitly using the orchestration boundary.

## Authorized native refactoring

The user explicitly authorizes rewriting SALMON source files for SALMONpy, not
merely wrapping the present CLI. Therefore communicator injection and explicit
state/context ownership take precedence over retaining global-variable layouts.
Refactor parallelization/communication and solver initialization at their proper
boundaries; avoid monkey-patching MPI_COMM_WORLD or repeatedly overwriting shared
globals as a substitute for independent state. Preserve numerical algorithms,
update order and the tested standalone entry point, not incidental file layout.

Apply full-solver changes in the externally prepared source and version them in
SALMONpy's pinned runtime patch/provenance. This authorization does not require
restoring the entire upstream source tree into this standalone package. The three
shared LD/SBE kernels retain byte equality. Embedded MPI support will be a new
capability with its own interface revision/build profile and MPI-compatibility
checks; current serial v2 bundles must not claim that capability.

## User-replaceable materials and A-to-J engines

This is a core requirement, not a later SALMON-specific extension. Built-in LD,
DFT/TDDFT and SBE adapters implement the same contracts available to user code.
The orchestrator must not branch on a closed list of material implementation
classes. Existing batch FDTD's closed LD/dielectric list remains legacy scope;
it must not be presented as the future plug-in API.

Two independently replaceable contracts:

1. `MaterialProvider.prepare(spec, context) -> InitialState` and
   `MaterialProvider.load(reference, context) -> InitialState`.
   Providers may run DFT, read a user tight-binding model, load experimental-fit
   parameters or construct another model. No universal SALMON input file,
   pseudopotential or wavefunction format is imposed. InitialState exposes an
   immutable descriptor and an opaque provider-owned payload. Descriptor fields
   include provider ID/version, state-format ID/version, units, symmetry/orientation
   if used, provenance, and the compatible response-backend family. Loading an
   incompatible state format is an error, never an implicit conversion.
2. `ResponseBackend.create_state(initial, point, step_config, context) -> State`
   and `State.advance(drive_frame) -> response_frame`, plus `close()` and declared
   capabilities. Storage, propagation algorithm and native implementation are
   private. Each call to create_state must create independent mutable state even
   when the immutable initial payload is shared. A memoryless A-to-J law is a
   valid special case; history-dependent responses retain their own state.

Provider and response factories are explicitly supplied Python objects or
registered entry points with package/version metadata. Do not evaluate strings
from input files to invent code, auto-download plug-ins or compile on import.
MaterialSpec config is provider-defined and schema-validated by the provider;
the core validates only common identifiers, units and coupling/resource metadata.
No inheritance from SALMON Fortran data types is required. Multiple providers and
response implementations may coexist in one MacroMesh.

Two execution adapters prevent unnecessary MPI complexity for users:

- `leader_only`: provider/response callback runs on the group's leader; framework
  performs status/result broadcasts and owns all group collectives. User receives
  no solver communicator and must not perform MPI collectives. Default group size
  is one; extra ranks are rejected unless a declared adapter has a useful role.
- `collective`: all group members call the backend in the same point order and
  receive the supplied solver communicator. Its communication and initialization
  must stay within that communicator. The backend explicitly declares whether a
  current is leader-complete, replicated or local contributions; the adapter
  produces exactly one group-reduced physical response. Avoid double reduction.

The public exchange contract always uses A/c and physical charge-current density
in declared atomic units, explicit times and point identity. Built-in SBE/TDDFT
raw-current conversion remains inside their adapters; user code already returning
physical current must not be negated again. Array dtype/shape, finite values,
step advancement and state ownership are enforced at the boundary. The user's
model is free to implement its own physics; SALMONpy cannot certify its physical
accuracy from interface compliance alone.

Capabilities include `execution_mode`, `multiple_states_per_group`, accepted
initial-state formats, supported coupling schemes, response time stages,
checkpoint support, and optional field-update contributions. A custom A-to-J
model does not have to implement SALMON's LD auxiliary Q or coefficient arrays.
It can be used with a compatible explicit/electronic coupling schedule. To claim
exact legacy semi-implicit LD coupling it must also supply the corresponding
optional prepare_field_update capability. Reject incompatible combinations at
planning time rather than silently changing the numerical scheme.

Proposed usage sketch (NOT executable/current public API):

```text
material = MaterialSpec(provider=user_provider, config=user_config)
response = user_response_factory
plan = plan_layers(layers, mesh, resources,
                   material_providers={...}, response_backends={...})
run(plan)  # creates separate state per point, supplies subgroup context
```

Add plug-in contract tests before exposing this API: custom scalar-memory response
against a serial reference; two points with distinct A histories from one initial
payload; built-in and custom materials mixed; swapped provider/material/point IDs;
leader-only versus collective test backend; wrong current sign/units declarations,
missing stage/point, exception before collective and unsupported capabilities.
Use manufactured laws as software tests, clearly separated from physical model
validation. The first working distributed example should include a small
user-defined Python response, so no hidden SALMON-only dependency can pass review.


## Implementation checkpoint — 2026-11-13

Completed: validated rank/point plans; solver/leader communicator splitting;
Fortran subgroup attach/detach in the real SALMON communication layer; optional
host-MPI bridge with source/object audits; independent user-created macro states;
point-keyed A broadcast and J gathering; LD endpoint-response adapter; callback
configuration/error/cleanup checks; paired guides and a shipped synthetic example.
Communicator owners reject premature close while native/session dependents live.

Validated locally: native groups of sizes 2 and 3; Python custom-response groups
of sizes 1/2/3 (6 ranks), two-material/four-point LD (4 ranks), and shipped example
(5 ranks). Full MPI SALMON rebuild preserves all 19 standalone FDTD interface
comparisons. Regression tests exercise mismatched plans, material-key collision,
noncontiguous user J, partial initialization cleanup and premature communicator
closure. The implemented fixed-endpoint response service is intentionally narrower
than the full DriveFrame/coupling-scheme protocol proposed above.

Still pending: explicit embedded DFT/TDDFT contexts and mutable-global extraction,
multiple electronic states per process, FDTD field-update/coupling phases,
layer meshing/resource estimation, persistence/restart and production scaling.
Do not infer these capabilities from the response-service demonstration.


## TDDFT service checkpoint

The experimental ResponseSession now has TDDFTMaterialProvider,
TDDFTPreparedState and TDDFTResponseBackend. This connects existing process-isolated
serial workers, one per macro point, and rejects multi-rank material groups.
Provider identity includes data/runtime/settings; the backend verifies serial
capability independently and launches the verified executable. Four point histories
and repeated runs match serial charge current exactly. Full embedded electronic
contexts/internal MPI remain pending; see 2026-11-13-tddft-context-inventory.md.

## Current allocation profile — user clarification

For the first electronic coupling implementation, one material and one macro
point per MPI solver group. Internal electronic MPI is still a future native
context milestone; the current process worker adapter uses one rank per group.
The standard electronic sample now uses two such groups and two points total.

For LD-dominated Maxwell calculations, use SALMON's already-parallel FDTD with
spatial domain decomposition. Each field rank owns many local LD cell/edge
histories and selects parameters by material ID; do not create a solver group
per LD cell. The Python FDTDMPI batch launcher exposes this existing capability.
It is distinct from the distributed prescribed-A response service and preserves
native LD implicit coefficients/Q ordering. See the local MPI validation record.
