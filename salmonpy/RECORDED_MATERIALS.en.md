# Connect recorded A-driven models to MPI

`bind_recorded_material(model, key=...)` reuses the existing endpoint-A adapter
while deriving its identity from `model.record()`. Supported examples are configured
SBE models and user Python models with explicit backend/new_state/timing contracts.
It does not run conformance tests or initialize MPI. Prepared models are borrowed:
close all response sessions before closing a model.

```python
from salmonpy.materials import bind_recorded_material, response_identity
from salmonpy.material_response import ResponseSession

binding = bind_recorded_material(model, key='prepared-film')
bindings = {10: binding, 90: binding}
point_map = {10: 'right', 90: 'left'}
record = response_identity(bindings, point_map, run_id='experiment-1')
session = ResponseSession(context, bindings, dt=0.02,
                          point_record=record, point_map=point_map)
```

All ranks use the same run ID/map, and must follow the collective lifetime rules
in [MPI_IDENTITY.en.md](MPI_IDENTITY.en.md). Integer transport order is sorted;
user IDs are preserved verbatim. Reuse one binding object for all points sharing
its cache key. Each point still creates its own mutable state. Each material group
must have exactly one MPI rank; this does not parallelize a serial kernel internally.

`record()` follows the [PointRun identity contract](POINT_RUNS.en.md): material ID
and definition SHA256, optionally both initial-state ID and digest. The binding
copies/freezes this record; direct MaterialBinding construction also validates and
freezes an optional material_record. Preparing and creating a new state rechecks
the model record. ResponseSession rejects disagreement between its point record
and a binding's material_record before callbacks. Definition hashes are declarations:
this does not detect a dishonest model that mutates its hidden state without changing
its record. Existing configured SBE loading checks its input files and owns the
loaded arrays; no background filesystem polling is added.

`response_identity` requires recorded bindings, validates the explicit integer map
and cross-point material/initial identity, and returns a detached point-run/1
record. It starts no workers and makes no MPI calls. The original bind_material
and unrecorded bindings remain supported by the original APIs.

LD remains E-driven; this helper does not silently convert A endpoints to E.
Use the explicit LDMaterialProvider/LDResponseBackend route where that particular
conversion is intended. TDDFT, including a ConfiguredMaterial wrapper around it,
is rejected here: use TDDFTMaterialProvider/TDDFTResponseBackend with the verified
managed runtime. A factory bridging configured TDDFT to that provider is still pending.
Collective native TDDFT and restart identity are outside this helper's scope.

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_recorded_material /new/output
mpiexec -n 3 python -m salmonpy.examples.mpi_recorded_material /new/sbe-output \
  --catalog /local/inputs/materials.json --material my-sbe
```

The first example uses a synthetic user model; the second selects a configured
SBE model and needs the installed native SBE wheel. Both compare MPI current with
the direct serial response and read saved identities back. No online data access
or compilation occurs. Exact equality is a regression observation, not a claim of
physical convergence or cross-platform bitwise equality.

Dedicated factories are now available separately: [CONFIGURED_NATIVE.en.md](CONFIGURED_NATIVE.en.md).
The generic bind_recorded_material still does not select TDDFT or convert LD implicitly.
