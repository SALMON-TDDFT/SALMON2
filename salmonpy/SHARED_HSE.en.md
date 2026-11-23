# Create independent fast-HSE calculations for the same material

[日本語](SHARED_HSE.md)

Use this when applying different light pulses to the same material, or placing
the same material at several macro points. Atomic structure and symmetry inputs
are shared; each calculation has its own wavefunctions, excitation history and
current history. Advancing one calculation does not change another.

Prepare the shared model with `NativeFastHSEBackend.open_model(initial, dt, context)`.
Use `model.new_state(...)` to start independent calculations from the original
ground state. Material-input checking/copying and MPI attachment occur once;
electronic initialization occurs on the first `new_state` call.

```python
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend

initial = FastHSEMaterialProvider(gs_directory, max_steps=100).prepare(context)
model = NativeFastHSEBackend(build_json, work_root).open_model(initial, .001, context)
states = []
try:
    with model:
        a = model.new_state('point-a')
        states.append(a)
        b = model.new_state('point-b', dt=.002, capacity=100)
        states.append(b)
        ja = a.step([0., 0., 1e-4])
        jb = b.step([1e-4, 0., 0.])
        # Use a distinct shared filename for each MPI rank.
        a.checkpoint(rank_checkpoint_path, False)
        restored = model.new_state('point-a')
        states.append(restored)
        restored.checkpoint(rank_checkpoint_path, True)
        j_next = restored.step([0., 0., 2e-4])
finally:
    for state in reversed(states):
        state.close()
    model.close()
```

`context` supplies `solver_comm` and `point_ids`; a caller-owned ParallelContext
can be used. All subgroup ranks call in the same order from the main thread.
The communicator must remain alive until every state/reference to a shared model is closed.
Input/library validation happens before the solver initialization is consumed; the
import itself loads no MPI or Fortran/C library and does not compile/download.

`new_state(point_id=None, *, dt=None, capacity=None)` defaults to the model's dt
and prepared maximum steps. Omitted point IDs get sequential integers. Explicit
string/integer IDs are retained on `state.point_id`; tokens are separate solver
identifiers. Each state owns A/current/clock/capacity, HSE histories and its own
output basename. `state.model_identity` identifies the common material snapshot.
Passing dt/capacity does not mutate the model or sibling states. The initial
source may be reused; keep the verified work copy unchanged.

The existing state methods `step(A)`, `advance(Aold,Anew)`, `checkpoint(path,False)`,
`checkpoint(path,True)` and `close()` retain their meanings. Direct load requires
a fresh state with matching dt/capacity and the same binary, GS/input and MPI
partition; it restores internal history without replaying accepted steps. The
calling program authenticates saved files/identity. Text logs start with the new state's
label; numeric restoration does not merge old text output. Ordinary protected
step journals select only that state's two streams, so rollback does not truncate
a sibling's output in the common directory.

`model.close()` prevents new creation; existing states can finish. It detaches
only after the last state closes. Closing all states while the model stays open
still permits fresh creation, without retaining a full hidden electronic state.
The model context manager does not close surviving states: always close them in
`finally`. A never-initialized closed model releases its unused Python reservation.
No MPI cleanup is triggered by Python garbage collection.

Scope of the legacy `open_model` route: one immutable **full-k or symmetry-reduced-k** model/process, existing k-only MPI1/2
calculation implementation and CPU/main-thread rules. Legacy `create_state` remains single-lifetime;
choose it or `open_model` for a process. An explicitly rebuilt runtime must expose
both handle protocol1 and model protocol1. The explicit heterogeneous registry below
adds different-input models. Automatic rank reassignment and threaded entry remain separate. With context_mpi_protocol1, initialization children retire at owner close.
Unscoped/model children remain attachment-owned until detach; repeated creation
is not an unlimited resource guarantee.
Budgets and physical-convergence advice remain advisory; inputs are not changed.
A is Cartesian A/c in atomic units, j is charge-current density (minus SALMON
matter current), dt atomic time and energy Hartree. x/y/z/circular drive is allowed.

Run [examples/shared_fast_hse.py](examples/shared_fast_hse.py) explicitly under the
matching MPI, with a prepared full-k zero-internal-pulse GS bundle:

```sh
mpiexec -n 2 python shared_fast_hse.py BUILD_JSON GS NEW_OUTPUT
```

The sample advances two polarizations and verifies fresh direct resume. It checks
composition/restart, not material convergence. See [solver handles](NATIVE_HANDLES.en.md)
for the underlying C ABI and remaining source ownership boundaries.

## Normal mixed MPI sessions

Opt in with `NativeFastHSEBackend(build_json, work_root, shared_model=True)`.
Pass `points_per_group={material_id: 2}` to `assemble_native_layers`; the
`ranks_per_point` value is the MPI width of each group. Every macropoint retains
its own state and ID. Points of the same material are grouped in encounter order,
including across separate layers. Omit the option to retain one point/group.
This is an explicit layout choice; budget warnings never change it automatically.
Known singleton TDDFT cannot share a group; user-written calculation implementations may provide their own
factory. Shared Fortran HSE calculation accepts full-k or an explicitly opted-in same-symmetry model.

`ResponseSession` optionally calls `backend.open_group(initial, dt, context)`
once per local binding. A returned factory supplies `new_state(point_id)` and
`close()`; `None` uses the existing `create_state` path. All subgroup ranks call
the same sequence. The session owns the factory, closes states before it, and
attempts cleanup of all created states even if one cleanup fails. Reusing a calculation implementation
object does not implicitly cache live sessions. User calculation implementations need no inheritance.

The installed normal sample supports `--hse-points-per-group 2`:

```sh
mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native \
  STANDARD_BUILD RUNTIME HSE_GS NEW_OUTPUT \
  --layers mixed_shared_hse_layers.json \
  --fast-hse-build SERIAL_HSE_BUILD --fast-hse-mpi-build MPI_HSE_BUILD \
  --hse-ranks 2 --hse-points-per-group 2 --steps 256 --dt .003 --protect-steps
```

Obtain the JSON with `importlib.resources.files('salmonpy.examples')`; use the
path to `mixed_shared_hse_layers.json`. It contains LD, SBE and two HSE points:
4 ranks with sharing, 6 without. `mixed_shared_tddft_hse_layers.json` adds PBE
(two internal ranks): 6 shared ranks or 8 separate ranks. Two field ranks need
one extra rank and `--field-grid 1 2 1`. Build explicitly and supply verified GS
bundles; the layouts do not download or manufacture a ground state.

Add `--direct-save-at 127`, then start a fresh job with
`--direct-resume SAVED_OUTPUT/direct-checkpoint` and the same layout/dt/steps.
Multi-binary groups use `salmonpy.native-direct/4` and per-rank, per-point files;
one-state filenames and schemas remain unchanged and versions1–3 remain readable.
Resume requires the same binary, inputs and MPI topology. Protected restore stages
independent candidates: a later-point failure preserves field, LD, SBE and all
electronic originals and allows retry without replay. This is not process-crash
recovery or an interchangeable cross-calculation implementation archive format.

Local Taylor and PT-CN tests compared shared versus separate MPI groups, independent
serial responses, and continuous versus direct resumed currents/fields/energies.
See the repository validation record for scope; this is not material convergence.

## Shared reduced-k symmetry

Use `NativeFastHSEBackend(..., allow_symmetry=True, shared_model=True)` in normal
sessions, or `allow_symmetry=True` with explicit `open_model`. A current runtime
must expose `sp_fast_hse_shared_symmetry_protocol()==1`; full-k remains compatible
with earlier model-protocol1 libraries. Every state uses the same immutable cell,
atoms, k mesh and retained operations. The density map is rebuilt synchronously
on construction; propagated wavefunctions, exchange maps/kernels and learned
histories are separate owners. Main-thread/collective ordering remains required.

For the normal mixed example add `--hse-symmetry --hse-reference-bundle FULL_GS`;
provide the reduced GS as the HSE bundle, and retain `--hse-points-per-group 2`.
The serial reference is full-k. For the explicit example add `--symmetry`; its two
drives are then along z (this still requires a group preserving z). Other drives
are allowed when the actual retained group preserves them; incompatible A is
rejected without modifying the accepted state. No projection or group conversion.
Taylor strict/learned and PT-CN off/learned were checked with internal MPI1/2,
independent controls, full/reduced-k and direct fresh restart. Different symmetry
groups within one model/process and cross-topology restart remain unsupported.

The supplied 3D localized mixed fixture was verified for128 steps (`--steps 128`,
`--direct-save-at 63`). At a longer256-step run its local field ceased to preserve
the retained group after145 steps, and input checking correctly rejected it while
protected stepping preserved the originals. A z-polarized source does not ensure
that every evolving local FDTD field stays along z. Use a symmetry-compatible
geometry/field or prepare full-k input for that calculation; do not discard small
components, rotate/project the drive or silently shorten a user's requested run.

A non-symmetric PBE point can generate components incompatible with the HSE
retained group even when the incident pulse is along z. The tested combined PBE
fixture is such a negative case; use full-k HSE for it. Full-k PBE/HSE mixtures
remain supported and were regression-tested with the new runtime.

## Per-context MPI initialization resources

A newly built runtime exposing `sp_fast_hse_context_mpi_protocol()==1` tags
communicators created during each synchronous state initialization with a unique
owner ID. `state.close()` releases its request caches, XC/exchange/FFT owners and
then those communicators. Other states and the calling program communicator stay live.
Registry slots are reused without reusing owner IDs during an attachment.
`sp_fast_hse_mpi_resources(int values[2])` reports live tracked children and registry
slots for explicit solver diagnostics (all subgroup ranks, main thread). It is
not an estimate of process memory or a general MPI resource count.

The shared density symmetry map uses MPI_COMM_SELF for the k-only, one-rank
r-space group so closing its last constructor does not invalidate siblings.
Other grid decompositions retain their supplied communicator. Unscoped children
remain attachment-owned; solver initialization crashes/allocation/MPI failures,
unrestricted configurations and concurrent entry remain outside this guarantee.
The explicit registry below additionally owns/swaps each model's density symmetry map.
Older runtimes still work in their documented scope but retain children until
detach. Rebuild explicitly for the new retirement behavior; never patch a binary.

## Different-input HSE models in one process

`backend.open_registry(context)` returns `NativeFastHSERegistry` on one existing
MPI managed by the calling program subgroup. `registry.open_model(initial, dt)` validates/copies each model
once; its first `new_state` registers input/unit/symmetry values. Subsequent states
start independently from that model's GS. The old `open_model` and `create_state`
routes remain unchanged. Choose a registry or one-initialization execution method for an attachment
lifetime; import never attaches/builds/downloads a solver runtime.

```python
backend = NativeFastHSEBackend(build_json, work_root, allow_symmetry=True)
with backend.open_registry(context) as registry:
    first = registry.open_model(initial_a, .001)
    second = registry.open_model(initial_b, .002)
    states = []
    try:
        a = first.new_state("material-a")
        states.append(a)
        b = second.new_state("material-b")
        states.append(b)
        ja = a.step([0., 0., 1e-4])
        jb = b.step([0., 0., -7e-5])
        a.checkpoint(rank_checkpoint_path, False)
        restored = first.new_state("restored-a")
        states.append(restored)
        restored.checkpoint(rank_checkpoint_path, True)
    finally:
        for state in reversed(states):
            state.close()
```

Inputs may have different cells, k representation, propagation/history
settings and symmetry tables. The existing Python provider still requires atomic
units (`unit_system=au`); owning unit tables does not lift that contract.
Each model owns its input values, unit tables,
complete symmetry workspace and density-symmetrization map. State selection restores
its model before step/checkpoint/transaction/close. Numerical/XC/exchange/FFT/MPI
owners remain per state; they are never copied as part of model selection.
All subgroup ranks select the same model and call on the main thread. The
process working directory is changed synchronously and restored by the adapter;
this is not thread-safe concurrent embedding.

Closing a model stops new states while existing states keep their model alive.
Closing the registry closes model factories; its solver attachment stays until
the last state closes. Close states explicitly, before calling program communicator/MPI
retirement. An idle model cannot release the registry's process reservation.
Collective model mismatch and work-directory failure are checked before solver
calls. General Fortran STOP/allocation/MPI/process failure recovery is unchanged.

Use a current explicitly built/audited runtime exposing
`sp_fast_hse_registry_protocol()==1`; old runtimes do not acquire this feature.
The public C registry ABI is in native_adapters/fast_hse_handles.h. C hosts supply
the selected model's working directory; positive state and negative model tokens
are opaque generation IDs. Same MPI subgroup/k-only MPI1/2 applies throughout;
changed-topology restore and automatic heterogeneous registry wiring in ordinary
material-ID/rank-plan assembly are separate. Existing same-model assembly stays.

Run the existing example with `--second-initial OTHER_GS --symmetry` to use two
prepared models. `--symmetry` permits the actual supplied symmetry declarations;
every retained group still checks the real drive, including x/y/circular A.
No projection/clamping is performed. Validation uses existing/synthetic GS inputs,
not converged material property benchmarks. Resource advice remains warning-only.

## Ordinary heterogeneous layer assembly

Use `NativeFastHSEBackend(..., registry_model=True)` in each HSE MaterialBinding
and `assemble_native_layers(..., model_groups={"hse-models": ["hse-a", "hse-b"]})`.
Each material retains its own provider/initial GS and every macrocell its own state.
Equal inner rank counts, one matching library/work root and matching symmetry
permission are required within the HSE group. Standard TDDFT and HSE still use
separate solver groups. Default plans and `shared_model=True` remain unchanged.

`salmonpy.examples.mixed_hse_models` demonstrates ordinary solver layer assembly,
FDTD-local LD, SBE, two HSE input IDs and protected direct save/resume. Its
`materials` JSON uses the existing `salmonpy.native-tddft-materials/1` local input
schema; this example explicitly chooses HSE providers, rather than guessing a
calculation implementation from filenames. Use k-only decomposition, e.g. `[2,1,1,1,1]`, per material.
The packaged `mixed_hse_models_layers.json` fixes layer IDs/order/thickness;
`--layers` accepts another layout. `--separate` creates independent material
groups for comparison. This40³/8³ small 3D fixture tests interfaces, not precision.

```sh
mpiexec -n 4 python -m salmonpy.examples.mixed_hse_models FIELD_BUILD HSE_BUILD MATERIALS_JSON NEW_OUTPUT
mpiexec -n 6 python -m salmonpy.examples.mixed_hse_models FIELD_BUILD HSE_BUILD MATERIALS_JSON NEW_OUTPUT --separate --compare GROUPED_OUTPUT
```

Use `--save-at 8` with total `--steps 32`; resume into a new output with
`--resume SAVED_OUTPUT/checkpoint --compare CONTINUOUS_OUTPUT`. Each material's
inner rank count and complete inputs must stay the same. Use a current runtime
with unique model output labels and collision-free symmetry input units.

The local-input template is `examples/mixed_hse_models_inputs.json`; replace bundle paths with your existing matching GS inputs.

## Per-call file resolution

With directory scope2/3, legacy/shared/registered fast-PR HSE calls use each work
folder as a solver file-resolution base without calling program chdir. Failure paths release
the scope. Older runtimes, including scope1, retain the prior route. This shares
GS call-input checking/cleanup code and does not enable concurrent threaded execution.
The numerical direct-checkpoint layout is unchanged.

Shared/registered fast-PR HSE states expose `state.timing_seconds` with optional timer protocol1: completed rank-local total, initialization, density and Hartree intervals in seconds. Query collectively on the material communicator. Each handle excludes parked, sibling and calling program work; unfinished intervals are absent, so total may be zero. Timing is not a numerical checkpoint payload, and restore does not erase spent work. Only the new query needs a matching runtime; older runtimes retain ordinary operations. Singleton timing has a separate legacy timer protocol1; threaded concurrency remains unsupported.

New runtimes also give the fixed-settings initialization guard an explicit model owner: one per registered model or same-material shared factory. Switching models no longer overwrites a process-wide guard. Only successful initialization commits it; failed creation preserves it. Validation predicates, C entry points and numerical checkpoint representation remain unchanged. The standalone legacy path retains its compatibility guard; this does not complete all shared-state isolation.

Reusing a closed model slot for a different material starts a fresh initialization guard for that slot, preserving every other live model.

Checkpoint array-size validation now uses call-local scratch, passed explicitly to each serializer. The existing state-derived bound, malformed-size checks and stream format remain unchanged. This does not add threaded concurrent solver execution.

Fast-PR HSE singleton runtimes with legacy timer protocol1 isolate initialization (with or without symmetry), propagation and close timing from the calling program. `state.timing_seconds` returns the same four completed counters as shared states. Old runtimes retain ordinary calculations; only the query requires this capability. Timing is not numerical checkpoint data, and MPI/update order/lifetime restrictions remain unchanged.

Checkpoint history g, b and coeff now use scalar IEEE inquiries over every element, retaining the existing metadata checks, rejection and format. Finite/NaN/positive- and negative-infinity tests use synthetic history; they are not new real-electronic-state disk-restart or Fujitsu compiler evidence.

`model_settings` is required. Every caller owns its settings, including the single-state and worker paths; there is no implicit process-global fallback. Process selectors still do not permit threaded reentry.

`model_settings` is required. Every caller owns its settings, including the single-state and worker paths; there is no implicit process-global fallback. Process selectors still do not permit threaded reentry.

## Fast-PR HSE performance mode

`NativeFastHSEBackend(..., execution_mode="performance")` selects the
caller-managed MPI route; the default is `"checked"`. The choice
is fixed on each state. Single-state, shared-model and heterogeneous-registry
routes are supported; serial subprocess workers are not supported by this option.

The performance update omits endpoint agreement, Python finite scans, electronic
and exchange-state backups, and numerical-recovery diagnostic votes/scans.
Initialization, library authentication, local lifetime/capacity/shape bounds,
explicit checkpoints and observations remain checked. Solver breakdown checks,
FFT/Hartree/exchange reductions and learning-readiness communication remain.
The mathematical compatibility of A with the retained symmetry little group
also remains: status34 rejects that input before mutation and permits retry.
Other returned failures mark the Python state failed, without an automatic
rollback guarantee; unhandled native failures can end the job.

Learning endpoint collection stays before energy/output handling in both modes.
An owner-local scheduling flag separates this ordering from recovery diagnostics;
it is reset after each update and is not an execution-mode selector in shared
process globals. Taylor/PT-CN, full/reduced k and internal1/2-rank synthetic
comparisons preserve current, energy and learning history, including explicit
in-memory restoration. These are not physical convergence or speed benchmarks.

`mpi_mixed_native --hse-execution-mode performance` requires an explicit
`--fast-hse-mpi-build`. It selects only HSE updates;
the field and outer coupling retain their independently selected checks.

## Modes of shared and registered fast-HSE states

A backend's `execution_mode` supplies its model default. Use
`model.new_state(point_id, execution_mode="checked")` or `"performance"`
to override per state; omit the keyword to inherit the model default.
`registry.open_model(..., execution_mode=...)` selects a registered model's
default. Session registry reuse preserves each material binding's chosen mode.
Siblings can mix checked/performance updates with separate electronic, exchange,
learning, A/current and output history. There is no process-global mode selector.

The performance handle update keeps local token generation/lifetime bounds and
the symmetry/A mathematical condition, while omitting diagnostic token votes
and automatic recovery copies. Physical activation of each model's input and
workspaces, solver communication and learning update order remain. Checked
siblings retain their full validation and recovery path. Explicit memory or disk
checkpoint operations retain checks/costs; restoring does not change the
receiving state's selected mode. Performance update failures still have no
automatic rollback guarantee.

`shared_fast_hse --execution-mode performance` demonstrates point-a in the
selected mode and point-b in checked mode, including a matching point-a restart.
With `--second-initial`, the same example uses two registered input models.
The usual `mpi_mixed_native --hse-execution-mode performance` also permits
shared HSE points; field/outer-coupling checks have independent settings.
