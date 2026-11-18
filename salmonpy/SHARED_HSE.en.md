# Python shared fast-PR HSE model

[日本語](SHARED_HSE.md) · Topic: `shared_fast_hse` · Revision: `2026-11-18.1`

`NativeFastHSEBackend.open_model(initial, dt, context)` returns a
`NativeFastHSEModel`. It verifies/copies one immutable full-k input bundle and
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

Scope: one immutable **non-symmetric full-k** model/process, existing k-only MPI1/2
backend and CPU/main-thread rules. Legacy `create_state` remains single-lifetime;
choose it or `open_model` for a process. An explicitly rebuilt runtime must expose
both handle protocol1 and model protocol1. Reduced-k/symmetry multi-context,
heterogeneous same-process inputs, automatic multipoint rank reassignment and
threaded entry are separate work. Child communicators remain attachment-owned
until final detach; repeated creation is not an unlimited resource guarantee.
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
