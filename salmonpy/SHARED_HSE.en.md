# Python shared fast-PR HSE model

[日本語](SHARED_HSE.md) · Topic: `shared_fast_hse` · Revision: `2026-11-18.5`

`NativeFastHSEBackend.open_model(initial, dt, context)` returns a
`NativeFastHSEModel`. It verifies/copies one immutable material/symmetry input bundle and
attaches one existing host-MPI communicator. The first `new_state` initializes
the native model; subsequent calls create fresh independent original-GS states.
There is no hidden numerical template and no per-state GS-directory copy.

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
The communicator must remain alive until every state/model lease is closed.
Input/library validation happens before the native lifetime is consumed; the
import itself loads no MPI or native library and does not compile/download.

`new_state(point_id=None, *, dt=None, capacity=None)` defaults to the model's dt
and prepared maximum steps. Omitted point IDs get sequential integers. Explicit
string/integer IDs are retained on `state.point_id`; tokens are separate native
identifiers. Each state owns A/current/clock/capacity, HSE histories and its own
output basename. `state.model_identity` identifies the common material snapshot.
Passing dt/capacity does not mutate the model or sibling states. The initial
source may be reused; keep the verified work copy unchanged.

The existing state methods `step(A)`, `advance(Aold,Anew)`, `checkpoint(path,False)`,
`checkpoint(path,True)` and `close()` retain their meanings. Direct load requires
a fresh state with matching dt/capacity and the same binary, GS/input and MPI
partition; it restores internal history without replaying accepted steps. The
host authenticates saved files/identity. Text logs start with the new state's
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
backend and CPU/main-thread rules. Legacy `create_state` remains single-lifetime;
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
composition/restart, not material convergence. See [native handles](NATIVE_HANDLES.en.md)
for the underlying C ABI and remaining source ownership boundaries.

## Normal mixed MPI sessions

Opt in with `NativeFastHSEBackend(build_json, work_root, shared_model=True)`.
Pass `points_per_group={material_id: 2}` to `assemble_native_layers`; the
`ranks_per_point` value is the MPI width of each group. Every macropoint retains
its own state and ID. Points of the same material are grouped in encounter order,
including across separate layers. Omit the option to retain one point/group.
This is an explicit layout choice; budget warnings never change it automatically.
Known singleton TDDFT cannot share a group; custom backends may provide their own
factory. Shared native HSE accepts full-k or an explicitly opted-in same-symmetry model.

`ResponseSession` optionally calls `backend.open_group(initial, dt, context)`
once per local binding. A returned factory supplies `new_state(point_id)` and
`close()`; `None` uses the existing `create_state` path. All subgroup ranks call
the same sequence. The session owns the factory, closes states before it, and
attempts cleanup of all created states even if one cleanup fails. Reusing a backend
object does not implicitly cache live sessions. User backends need no inheritance.

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
recovery or an interchangeable cross-backend archive format.

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
the retained group after145 steps, and admission correctly rejected it while
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
then those communicators. Other states and the host communicator stay live.
Registry slots are reused without reusing owner IDs during an attachment.
`sp_fast_hse_mpi_resources(int values[2])` reports live tracked children and registry
slots for explicit native diagnostics (all subgroup ranks, main thread). It is
not an estimate of process memory or a general MPI resource count.

The shared density symmetry map uses MPI_COMM_SELF for the k-only, one-rank
r-space group so closing its last constructor does not invalidate siblings.
Other grid decompositions retain their supplied communicator. Unscoped children
remain attachment-owned; native initialization crashes/allocation/MPI failures,
unrestricted profiles and concurrent entry remain outside this guarantee.
The explicit registry below additionally owns/swaps each model's density symmetry map.
Older runtimes still work in their documented scope but retain children until
detach. Rebuild explicitly for the new retirement behavior; never patch a binary.

## Different-input HSE models in one process

`backend.open_registry(context)` returns `NativeFastHSERegistry` on one existing
host-MPI subgroup. `registry.open_model(initial, dt)` validates/copies each model
once; its first `new_state` registers input/unit/symmetry values. Subsequent states
start independently from that model's GS. The old `open_model` and `create_state`
routes remain unchanged. Choose a registry or legacy route for an attachment
lifetime; import never attaches/builds/downloads a native runtime.

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
Closing the registry closes model factories; its native attachment stays until
the last state closes. Close states explicitly, before host communicator/MPI
retirement. An idle model cannot release the registry's process reservation.
Collective model mismatch and work-directory failure are checked before native
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
