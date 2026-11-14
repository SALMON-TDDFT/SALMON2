# Author one material model for serial and MPI use

[日本語](MATERIALS.md) · Topic: `materials` · Paired revision: `2026-11-14.2`

Implement one serial material contract. Use `crystal.respond` for a waveform,
and `materials.bind_material` to connect it to the existing MPI response service.
MPI distributes material points, each with independent mutable history. This
adapter does not parallelize one serial state across several ranks.

## Model contract

The initial model is read-only and exposes a nonempty algorithm/version `backend`.
`new_state(*, dt, time0=0., initial_A=(0,0,0))` returns fresh independent state.
Do not mutate shared initial arrays or incoming A arrays.

| State entry | Contract |
|---|---|
| `time` | Accepted time, initially time0, then time0+n·dt |
| `current` | Finite three-component current density at accepted time |
| `timing` | Explicit TimeSpec: fixed dt/time origin, endpoint A and current |
| `step(A_new)` | Accept next endpoint A and return current there; update current consistently |
| `close()` | Release resources; optional for states without resources |

A means SALMON Ac=A/c. Return laboratory-frame **charge-current density j** in
atomic units. Initial current must be physical rather than a placeholder.
Declare `TimeSpec(dt, (InputPort('A'),), (OutputPort('j'),), t0=time0)`.
If the internal observable is SALMON matter current J, explicitly declare
`OutputPort('J_matter')`. Both serial respond and this adapter then convert
j=−J exactly once. The software cannot infer physical units/signs from numbers.

Check calculations before committing state so a failing step leaves accepted
history intact. A model calling MPI internally must use the collective backend
contract below. Models/states are not safe for concurrent calls from threads.

## Use the same model serially

[custom_material.py](examples/custom_material.py) contains an executable example.
Replace its model/state classes with your own constitutive law or library calls.
Its memory recurrence is a manufactured interface test, not a real-material fit.

```python
from salmonpy.examples.custom_material import ExampleMaterial, waveform
from salmonpy.materials import check_material, bind_material
from salmonpy.crystal import respond

model = ExampleMaterial(gain=2., decay=0.8)
drive = waveform()
report = check_material(model, drive)
response = respond(model, drive)
binding = bind_material(model, key='sample', identity=model.identity)
```

`check_material(initial, drive, rtol=1e-12, atol=1e-14)` uses only the caller's
explicit representative A waveform. It creates two states and checks finite
values, timing, input ownership, independence and repeatability. E samples and
clearly nonuniform grids are rejected before creation; cumulative clock consistency
is also checked during evolution. Success returns immutable
MaterialCheck. Its states are closed on success/failure; the caller's model is
never closed. Start with a short drive that fits two states and their histories.

It does not prove physics, units, all possible inputs, complete model immutability,
failure atomicity or MPI collectives. These appear in report.not_tested.
rtol/atol compare repeated identical algorithms, not physical accuracy.
Binding or starting a production session never implicitly runs this checker.

## MPI binding

The result is an existing MaterialBinding passed directly to ResponseSession.
Reuse the same binding instance for points with one key. The author must supply
identity identifying the algorithm version, all parameters and initial data.
Do not use Python addresses or rank numbers. Incorrectly assigning the same
identity to different settings may hide differences; this is a declaration, not
an automatic inspection of arbitrary Python objects. Keep model configuration
and identities unchanged during a session.

```sh
python -m salmonpy.examples.custom_material
mpiexec -n 2 python -m salmonpy.examples.custom_material --mpi NEW_JSON
mpiexec -n 3 python -m salmonpy.examples.custom_material --mpi ANOTHER_JSON
```

The MPI example has one field rank and one rank per material group, distributing
four points. A group may own several points. Initial A is zero, matching the
existing service. Each complete point current is collected once, without rank
or point multiplicity. Importing the module neither imports/initializes MPI nor
launches ranks. Only MPI execution needs the optional salmonpy[mpi] environment
matched to the MPI launcher. See [DISTRIBUTED.en.md](DISTRIBUTED.en.md).

The new adapter requires explicit timing, without guessing missing declarations
from legacy models. check_material uses waveform time0/initial A; the MPI service
starts with time0=0 and A=0. Implement the conditions needed for both uses.

Contract errors during creation/stepping raise. A failed adapter cannot be reused.
An MPI step failure invalidates the session collectively; other points may already
have advanced, so no global rollback is promised. Collectively close the session
before its context.

## Internally parallel and native models

A genuinely parallel model continues to implement the existing provider/backend
collective contract using solver_comm. This serial adapter rejects multirank
material groups and does not duplicate work while claiming internal parallelism.
See [DISTRIBUTED.en.md](DISTRIBUTED.en.md).

Native TDDFT retains TDDFTMaterialProvider/TDDFTResponseBackend with managed-runtime
validation. LD remains E-driven and uses LDMaterialProvider/LDResponseBackend.
This entry does not implement whole-MPI rollback, collective checkpoint restore
or an MPI-wide C ABI. For serial native material plugins, see [C_ABI.en.md](C_ABI.en.md).
