# Checkpoint capability and compatibility

[日本語](CHECKPOINTING.md) · Topic: `checkpoint_compatibility` · Paired revision: `2026-11-17.1`

`salmonpy.checkpointing` provides a common admission step for built-in and user
backends. Call it before your own SAVE/LOAD functions. It returns detached JSON
metadata after checking all participating ranks. Serial hosts pass `comm=None`;
MPI hosts supply their existing communicator. No solver context or inheritance
is required. Import and checking perform no MPI initialization, file access,
compilation, download, propagation or persistence callbacks.

## Declare capability

```python
from salmonpy.checkpointing import CheckpointCapability

capability = CheckpointCapability(
    "my-material-json/1", ("save", "restore"), (1,),
    identity="implementation-and-fixed-model-input-digest",
)
record = capability.to_dict()
assert CheckpointCapability.from_dict(record) == capability
```

The object is immutable. `format` identifies a versioned payload format;
`modes` is a tuple containing save and/or restore. `decomposition` is a positive
integer tuple whose product equals the actual group size. Its axis order is part
of the contract: `(2,1)` differs from `(1,2)`. Pin implementation, fixed model/input
and runtime identity with `identity`. An empty identity is omitted from metadata
for compatibility with existing native-direct/2 archives. Checks compare the
explicit declarations; model identity and file integrity remain caller-owned.

## Serial user backend

The packaged example contains a complete user model, JSON persistence, payload
SHA256 verification, independent states and comparison with continuous execution:

```sh
python -m salmonpy.examples.custom_checkpoint NEW_OUTPUT
```

It saves at17 and advances only15 more steps through32; current is exactly equal
to continuous execution. Changing the model is rejected before its loader runs.
The recurrence is synthetic interface data, not a physical material prediction.

This is the admission part of that workflow:

```python
from salmonpy.checkpointing import check_checkpoint_compatibility
from salmonpy.examples.custom_checkpoint import UserCurrent

groups = [{"name": "user", "role": "material", "ranks": [0], "point_ids": ["film"]}]
state = UserCurrent()
saved = check_checkpoint_compatibility(
    None, groups, states={"film": state.checkpoint_capability}, operation="save")
# The host now saves its payload and the returned metadata, and verifies hashes.
candidate = UserCurrent()
check_checkpoint_compatibility(
    None, groups, states={"film": candidate.checkpoint_capability},
    operation="restore", saved=saved,
    required_files=("state.json",), available_files=("state.json",))
# After successful admission, the host can call candidate.load(its_verified_path).
```

`available_files` is a names-only inventory. The API checks safe relative POSIX
paths and required-file presence; it does not read files or verify their hashes.
The host must authenticate payloads and check that the restore candidate is fresh
before calling its own loader. The sample above illustrates metadata admission;
the runnable example supplies the actual I/O and checksum checks.

## MPI hosts and failure agreement

All ranks of the supplied communicator call the API in the same sequence.
`groups` is a list/tuple of ordinary mappings containing name, ranks and optional
role/point_ids. Ranks are relative to this communicator; every rank belongs to
exactly one group. Material groups own distinct integer or string point IDs
(including distinct string representations). Field groups use role=field and no
point IDs. Standalone material groups need no field group.

Each rank supplies only its own group's `states` mapping and, on field ranks,
its field capability via `field`. Capabilities can be objects or explicit mapping
records. Group size, point coverage, within-group declarations, complete rank
ownership and the saved compatibility record are checked together. All ranks
must agree on operation and saved metadata, including unchanged group rank order.

Use plain declaration mappings for per-rank unvalidated configuration so errors
are collected inside the call. If constructing/inspecting a capability fails
before the call, catch the exception and pass its text as `local_error` while
still entering the collective. An error on one rank then prevents every rank
from proceeding to LOAD. The caller must supply a live communicator; it retains
ownership of initialization, splitting and finalization. The API accepts its
rank/size/allgather interface and imports no mpi4py itself.

```python
from mpi4py import MPI
from salmonpy.checkpointing import check_checkpoint_compatibility

comm = MPI.COMM_WORLD  # Example requires three ranks, already initialized by host.
groups = [
    {"name": "field", "role": "field", "ranks": [0], "point_ids": []},
    {"name": "material", "ranks": [1, 2], "point_ids": ["film"]},
]
field = {"format": "user-field/1", "modes": ["save", "restore"], "decomposition": [1]}
material = {"format": "user-current/1", "modes": ["save", "restore"],
            "decomposition": [2, 1], "identity": "fixed-model-and-runtime-digest"}
saved = check_checkpoint_compatibility(
    comm, groups, operation="save", field=field if comm.rank == 0 else None,
    states={} if comm.rank == 0 else {"film": material})
```

For restore, pass operation=restore and the same saved metadata on all ranks.
Required-file inventories may be rank-local. Errors use `salmonpy.errors.InputError`
(a ValueError), and successful calls return `salmonpy.direct-compatibility/1` JSON.

## Existing mixed native calculations

The native mixed sample now uses this common API, retaining solver-specific
fresh-state and native-ABI checks. Existing native-direct/2 compatibility metadata
is unchanged; pinned native payloads and numerical update order are unchanged.
Custom direct-state/1 sample adapters may expose `checkpoint_capability` as an
object/mapping. The previous checkpoint_format/checkpoint_modes/
checkpoint_decomposition attributes remain accepted. The sample's callback/clock/
A/current requirements are described in [COUPLED.en.md](COUPLED.en.md).

Current admission requires the same topology and declared format/identity. Payload
conversion, changed-topology redistribution, generic archive I/O and atomic rollback
remain backend/host responsibilities. The API does not execute or certify an
arbitrary backend's persistence implementation. Historical native-direct/1 archives
without compatibility metadata require their original matching reader.
