# Checkpoint capability and compatibility

[日本語](CHECKPOINTING.md) · Topic: `checkpoint_compatibility` · Paired revision: `2026-11-17.1`

`salmonpy.checkpointing` provides common admission and execution for built-in and user
backends. `execute_checkpoint` calls explicitly bound SAVE/LOAD functions after admission. It returns detached JSON
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

## Execute backend-owned persistence

Bind zero-argument callbacks with `CheckpointTarget`, then call `execute_checkpoint`
with operation=save or restore. Paths and payload formats belong to the backend;
use closures or `functools.partial` to supply arguments. A complete file-based
example runs with `python -m salmonpy.examples.custom_checkpoint NEW_OUTPUT`.

The same entry can also transport an in-memory backend token:

```python
from salmonpy.checkpointing import CheckpointCapability, CheckpointTarget, execute_checkpoint

groups = [{"name": "user", "ranks": [0], "point_ids": ["film"]}]
cap = CheckpointCapability("user-memory/1", ("save", "restore"), (1,), identity="model-v1")
state = {"value": 17}
saved_run = execute_checkpoint(
    None, groups, operation="save",
    states={"film": CheckpointTarget(cap, save=lambda: dict(state))})
payload = saved_run["states"]["film"]
candidate = {}
execute_checkpoint(
    None, groups, operation="restore", saved=saved_run["compatibility"],
    states={"film": CheckpointTarget(cap, restore=lambda: candidate.update(payload))})
assert candidate == state and candidate is not state
```

The result contains `compatibility` JSON, a process-local `field` callback result
and a `states` mapping keyed by string point IDs. Opaque tokens, arrays and receipts
remain local; they are not copied or sent through MPI. The backend controls their
ownership. Callback return values, including None, are accepted; callbacks must
raise exceptions to report failure. A return value of False is not an error signal.

Capability and callback availability are checked collectively before the first
callback. Groups execute their own point order in synchronized slots, so groups
with different point counts are supported. An ordinary callback exception becomes
RuntimeFailure on every rank; later slots stop. Completed state changes remain.
A collective backend must complete its own MPI protocol before raising; the common
entry cannot unblock ranks stranded inside a backend collective or recover a
Fortran STOP/process crash. Hosts decide whether to discard affected candidates
or abort their native job. MPI initialization, splitting and finalization stay
with the host. Metadata/hash authentication and completion publication also remain
host/backend policies; the native sample preserves its existing inventory/marker.

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
before binding its own loader to the execution entry. The sample above illustrates metadata admission;
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

## Preserve live state on restore failure

Use `restore_checkpoint` for protected restoration. Keep every mutable local
component that must switch together in one `CheckpointState` holder. Pause
computation during the call, and access the states through the holder afterward.
`stage_restore()` loads and validates a fully independent candidate; it must never
mutate original states or share mutable numerical resources with them. Returning
the original object, a reused candidate or None is rejected.

```python
from salmonpy.checkpointing import (CheckpointCapability, CheckpointTarget,
                                   CheckpointState, restore_checkpoint)

groups = [{"name": "user", "ranks": [0], "point_ids": ["film"]}]
cap = CheckpointCapability("user-memory/1", ("save", "restore"), (1,), identity="model-v1")
from salmonpy.checkpointing import check_checkpoint_compatibility
saved = check_checkpoint_compatibility(None, groups, states={"film": cap})
original = {"step": 5}
live = CheckpointState(states={"film": original})
run = restore_checkpoint(
    None, groups, live=live, saved=saved,
    states={"film": CheckpointTarget(cap, stage_restore=lambda: {"step": 17})})
assert live.states["film"]["step"] == 17
assert run["previous"]["states"]["film"] is original
assert original == {"step": 5}
```

The common entry stages all components, shares errors across MPI, then votes before
publishing each rank's entire local bundle with one pointer assignment. Publication
calls no user code. Previous states remain strongly referenced in `run['previous']`;
retire their resources explicitly after success. Existing external aliases continue
to refer to old objects, so recreate bindings or use closures reading the holder.
The operation uses additional memory/resources while old and candidate states coexist.

On an ordinary staging failure, `AtomicRestoreError` (a RuntimeFailure) reports
`original_preserved=True`. All holders retain the original bundle, and computation
can continue from the originals. Optional `discard_restore(candidate)` cleans
attempted candidates in reverse synchronized slots. It is called on every relevant
rank, including `discard_restore(None)` when staging threw or returned an unsafe
alias. It must handle partial staging resources and complete its subgroup protocol
before raising. A stage factory can instead clean its own failed allocations before
raising. Cleanup errors are reported in `error.cleanup_errors`, separately from
preservation of the live states. Candidate cleanup never authorizes changing live
objects. Under the backend isolation contract, no original is passed to discard.

`execute_checkpoint(..., operation='restore', atomic=True, live=live)` is the same
protected path. Existing execute_checkpoint defaults to non-atomic execution for
compatibility; direct-only targets are rejected by the protected path with zero
LOAD calls, without fallback. LD/SBE states can be staged via new_state plus
snapshot restore; custom JSON factories are shown in the runnable example.

### Native numerical candidates (experimental)

Explicitly rebuilt native FDTD, standard TDDFT and fast-PR HSE expose
checkpoint_protocol=1: modes0 save,1 legacy immediate restore,2 stage a detached
numerical candidate,3 publish prepared memory,4 discard. One candidate per native
lifetime is allowed. Pause evolution, finish mode2 and agree success across every
host group before mode3. On handled failure every group discards with mode4 and
can continue its original. HSE retains validated exchange/history data in memory;
mode3 never reopens the checkpoint. Reject older libraries via the protocol query.

The internal `salmonpy._native_checkpoint.NativeRestoreCandidate` helper and
`tools/test_native_restore_candidates_mpi.py` demonstrate the common entry:
bind stage/discard callbacks, use the common API for collective preparation,
then call publish in each group and install preprepared Python A/current/clocks.
The generic holder switch alone does not finish native publication. This is an
internal experimental contract, not a stable public C ABI.

Only mutable numerical candidates are detached. Initialized MPI/Libxc/FFT resources
remain in the existing context; this does not enable concurrent native lifetimes
or continued use of an old alias after publication. Publication uses prepared
memory but includes allocating assignments: allocation failure/process death is
outside the recoverable contract. The ordinary mpi_mixed_native integration is described below.
Exact-layout/input/hash admission is still required.

This contract covers handled exceptions with healthy participating MPI ranks and
correct isolated factories. It does not recover process crashes, communicator loss,
Fortran STOP or a backend which mutates/shares its original state during staging.

## Existing mixed native calculations

The native mixed sample uses both common admission and execution, retaining solver-specific
fresh-state and native-ABI checks. Existing native-direct/2 compatibility metadata
is unchanged; pinned native payloads and numerical update order are unchanged.
Custom direct-state/1 sample adapters may expose `checkpoint_capability` as an
object/mapping. The previous checkpoint_format/checkpoint_modes/
checkpoint_decomposition attributes remain accepted. The sample's callback/clock/
A/current requirements are described in [COUPLED.en.md](COUPLED.en.md).

Current admission requires the same topology and declared format/identity. Payload
conversion, changed-topology redistribution, generic archive I/O and atomic rollback
remain backend/host responsibilities. The execution entry runs supplied functions without certifying the backend's
persistence implementation or atomic rollback. Historical native-direct/1 archives
without compatibility metadata require their original matching reader.

## Protected restore in the ordinary mixed sample

`mpi_mixed_native --direct-resume PATH` now defaults to protected restore.
Initialized field/standard TDDFT/fast-PR HSE require candidate-capable runtimes.
The field candidate includes PML and every LD material's P/J/update history.
SBE loads snapshots into independent new handles. Python clocks, currents, A,
response adapters, ownership lists and frames.json are prepared before publication.
All MPI participants must succeed before anything live is published. No extra SCF
or time propagation is performed by restore.

The internal `examples._direct_checkpoint.checkpoint(..., restore=True)` also
accepts advanced paused sessions. A staging failure raises AtomicRestoreError on
all ranks and preserves original field/LD/electronic states, clocks, references
and failed flags. Make the same collective decision everywhere: continue the old
trajectory or remove the failure cause and retry the same checkpoint. Input/history
admission errors are rejected before LOAD as well.

After success retrieve current states through session.states; old SBE handles are
retired. Release errors are separate coupled.checkpoint_cleanup_errors diagnostics,
not a failed restoration. Native resources still have one lifetime; allocation
failure during publication, process death and failed MPI remain outside recovery.

Older ABIs, process-isolated HSE workers and custom direct-only LOAD targets are
rejected before loading. Only explicitly request legacy behavior with CLI
`--unprotected-direct-restore` or internal `protected=False`; it does not preserve
original states on mid-load failure. There is no automatic fallback.

`tools/test_mixed_protected_restore_mpi.py` takes the ordinary sample arguments,
injects a rank failure after loading a later electronic candidate, then checks
original preservation (including two field-local LD materials), continued old
evolution and exact retry results for SBE plus internal-MPI TDDFT/fast-PR HSE.
Source-tree tests can explicitly select a built LD/SBE library with `--kernel LIB`;
installed-wheel users do not need that override.

Live rewind does not truncate native output files already written. Manage the
trajectory with the saved frames plus subsequent step returns, or resume in a
fresh MPI job/output directory to keep native output streams separate.

Memory advice includes concurrent candidate storage; see [MEMORY.en.md](MEMORY.en.md). Exceeding the ordinary sample budget only warns and does not refuse protected restore.

## Mixed sample identity and restart

`python -m salmonpy.examples.mpi_mixed_native` preserves `--run-id`, logical point IDs, material IDs, material-definition and initial-state SHA256 digests, MPI ownership, and the LD cell-material mapping. The default `mixed-example` is a fixture label; provide an explicit run ID for each experiment. Output and temporary directories are excluded from definition fingerprints.

The same identity appears in `identity.json`, normal-run `result.json`, and the `point-response` archive. Identified direct checkpoints use `salmonpy.native-direct/3`; replay checkpoints use `salmonpy.native-replay/2`. Numerical payload formats and update order are unchanged. Direct restore checks identity and MPI ownership before loading numerical state; save checks agreement across all ranks.

Use `--legacy-checkpoint-identity` explicitly only for old unidentified checkpoints. This cannot bypass identity mismatches in identified checkpoints and does not disable protected restore. Custom hosts supply `ResponseSession` point records through `point_record` / `point_map`. These are caller declarations, not certification of physical accuracy.

## In-memory snapshots with point identity

`PointState.snapshot()` wraps the backend-owned snapshot with immutable `salmonpy.point-snapshot/1` identity. `state.restore(saved)` checks run, point, material-definition and initial-state identity plus package version before calling the backend. Cross-point/run and raw backend tokens are rejected. Even with explicit matching run/point IDs, the backend must still validate its model, time step and other settings.

This is process-local rewind, not branching, cloning, disk persistence or an MPI-wide transaction. Unsupported custom backends raise `SnapshotError`. Backend restore-failure guarantees remain the backend's responsibility; identity rejection invokes no backend restore callback. Raw-state snapshot APIs are unchanged.

`python -m salmonpy.examples.point_ld /new/output/directory` also verifies cross-point rejection and exact current continuation after same-point restore.

## Native TDDFT rejection before update (limited scope)

Newly prepared and built standard embedded TDDFT exposes `sp_tddft_error_protocol() == 1`. Validation is collective on the attached material communicator:10 invalid/nonfinite input, stale A or missing input file;11 rank disagreement;12 unavailable lifecycle or repeated initialization;13 exhausted step capacity;14 pending restore candidate. These codes return unanimously before numerical mutation. Initial admission codes10/11 do not consume the lifetime; Python initialization can be retried with corrected input. Failed work directories remain for diagnosis.

Python maps only protocol-declared unanimous rejection to `InputError`. Old libraries, unknown codes or nonfinite results retain `RuntimeFailure`, with no inferred state guarantee. Every solver-group member must enter the same call. Invalid pointers, MPI/process failure, allocation errors, deep solver STOP and failures after numerical updates are outside this admission-only contract; see the limited deep-recovery contract below. Actual initialization still consumes the one-process lifetime. Fast HSE and native FDTD error contracts are unchanged.

This guarantee concerns an individual TDDFT state's pre-update rejection. A coupled session may already have advanced another component; it remains failed under the existing driver policy. This does not implement general coupled-step rollback.

## Direct rollback after deep TDDFT numerical errors

New standard native TDDFT exposes `sp_tddft_recovery_protocol() == 1` and retains a per-step numerical payload in memory. Collective checks cover propagated wavefunctions(20), density(21), Hartree field(22), XC potential/energy(23), current/total energy(24), and a protected propagator error(25). Every material-group rank agrees before returning. Errors restore accepted orbitals, densities, fields, propagation history, XC state, counters and original overlap flags without propagation replay. Failed attempts add no ordinary RT output records.

A state returned by Python `NativeTDDFTBackend` raises `NumericalError` from `advance(A_old,A_new)` only for protocol-declared unanimous recovery. Its `A` and `steps` remain accepted; correct the cause and retry from that time. No automatic pulse weakening or dt changes are applied. Old libraries and unknown errors retain conservative `RuntimeFailure` handling without inferred recovery.

Supported: periodic, unpolarized, fixed-ion middlepoint PZ, automatic/fixed TB-mBJ, Libxc-PZ/PBE. SOC,+U,jellium, ionic motion/force output, Gram–Schmidt, projection/unfold and microscopic Maxwell coupling are outside the protected profile and return15 before updating. Ordinary batch calculations retain their existing route; external checkpoint encoding is unchanged. Per-step copies/checks add memory and runtime overhead: leave resource headroom. Budgets remain warning-only.

This is a TDDFT-group guarantee, not automatic rollback of an entire mixed FDTD step whose other components may have advanced. MPI/process death, allocation failures, uncontrolled third-party STOP, deep initialization STOP and output-device failures are excluded. One-rank fault injection exists only in a separately compiled `SALMONPY_RECOVERY_TESTING` build; production builds export no injection entry.
