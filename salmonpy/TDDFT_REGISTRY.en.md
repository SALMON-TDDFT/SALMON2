# Different standard TDDFT models in one process

[日本語](TDDFT_REGISTRY.md) · Topic: `heterogeneous_tddft_models` · Revision: `2026-11-18.1`

`NativeTDDFTBackend.open_registry(context)` creates one explicit host-MPI
attachment. `registry.open_model(initial, dt, decomposition=...)` registers an
independent immutable material model; `model.new_state(point_id)` returns the same
A/c→charge-current state interface as the shared-model API. Different models can
have different inputs, XC, cells/grids/k meshes/atoms/PP references, dt/capacity and
inner decomposition. All use the same communicator; each decomposition product
must match its size. A model's own inputs are fixed, and each state starts its
original GS independently. Material/point labels belong to the host, not opaque
native tokens.

```python
from salmonpy.embedded_tddft import NativeTDDFTBackend

# context.solver_comm is the host's initialized MPI subgroup. Prepare initial_a/b
# with NativeTDDFTMaterialProvider and verified local GS/PP/input bundles.
backend = NativeTDDFTBackend(build_record, work_root,
                            decomposition=(1, 2, 1, 1, 1))
registry = backend.open_registry(context)
models, states = [], []
try:
    a = registry.open_model(initial_a, 0.003)
    models.append(a)
    b = registry.open_model(initial_b, 0.004, decomposition=(1, 1, 1, 1, 2))
    models.append(b)
    sa = a.new_state("film-a/0")
    states.append(sa)
    sb = b.new_state("film-b/0")
    states.append(sb)
    ja = sa.step([1e-4, 0.0, 0.0])
    jb = sb.step([0.0, 1e-4, 0.0])
    # Calls can alternate; sa.time uses a.dt, sb.time uses b.dt.
finally:
    for state in reversed(states):
        state.close()
    for model in reversed(models):
        model.close()
    registry.close()
```

Choose decompositions compatible with each saved state. Differing dt demonstrates
independent clocks; it does not add FDTD subcycling. The ordinary layer sample
still assigns different material IDs to separate groups. Automatic placement of
different material bindings in one shared response group is separate integration
work. Explicit registry hosts can use the states directly.

Use an explicitly built/verified CPU runtime with `sp_tddft_registry_protocol()==1`.
The profile remains periodic, full-k/no symmetry reduction, unpolarized fixed-ion
PZ/TB-mBJ automatic/fixed-c/Libxc-PZ/PBE, middlepoint transverse, initial A=0 in Python.
x/y/z and rotating A work. A is Cartesian A/c, dt atomic time, output j is atomic
charge-current density; do not negate it again. HSE is the separate fast-PR engine
and is not registered here. Imports do not initialize MPI, build or download.
Keep model bundles/work files unchanged; one work copy per model is made in the
shared writable work root, not per state. Registering initializes a first state
that is held until first new_state (or closed if unused).

All subgroup ranks call in the same order on the main thread. Registry mode is
explicit and cannot mix legacy initializers or implicit model_new. Default
single-model entry points remain unchanged. Model/registry close prevents new
creation while live states can finish; close all states before host MPI finalizes.
Closing one model leaves siblings intact. Group-wide model choices are checked
before queued state delivery, and model cleanup uses deterministic token order.
A handled close failure attempts remaining owners and permits close retry.

Each model owns parsed/default inputs, a post-initialization value view, initial A,
optimization modulo tables/flags and XC scratch. Array ownership is parked/switched
with move_alloc, without copying wavefunctions, Libxc handles or MPI owners. Each
state owns wavefunctions/density/XC payload, clock/A/current, streams and staged/
protected state. PP format flags reset during species initialization. Timers,
attachment and synchronous compatibility views remain process infrastructure;
this is sequential switching, not concurrent threads or unrestricted solver reentry.
General parser/allocation/I/O/process/MPI failure recovery is not promised.

`state.checkpoint(path)` and fresh-state `checkpoint(path, restore=True)` work per
model, without history replay; filenames are rank-local. Authenticate the same
binary/inputs/model/dt/capacity/decomposition and files before loading: raw files
are not a portable archive. Existing owner-specific protected-step/restore helpers
use the same positive state handles. Automatic ordinary-plan aggregate registry
construction and changed-topology portable restart are separate work.

Example: `python -m salmonpy.examples.tddft_models BUILD INITIAL_A INITIAL_B NEW_OUTPUT`
under matching MPI ranks (default orbital split). Developer validation tool
`tools/test_heterogeneous_tddft_mpi.py NEW_BUILD OLD_BUILD INITIAL_A INITIAL_B OUTPUT
--size 2 --b-split rgrid` uses six ranks: two models in the first group and two
independent old-runtime controls in the others. It checks interleaved trajectories,
direct resume, closing siblings, idle recreation, rejected/stale/mixed tokens,
failed registration and rank-independent cleanup. Existing Si fixtures with
16³/8³ grids, different cells/atom/k counts and five XC choices check ownership;
no new GS is computed. This is implementation parity, not material convergence.
See [shared TDDFT](SHARED_TDDFT.en.md) for base states and parser ownership, and
[native registry header](native_adapters/standard_tddft_models.h) for C/C++.
