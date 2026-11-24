# Checkpoint capability and compatibility

## Start from what you want to save

| Existing calculation | Read next |
|---|---|
| LD, SBE or a field state; temporary in-memory rewind | [snapshot / restore](SNAPSHOTS.en.md) |
| Retained GS; resume SCF with mixing history | [SCF save/resume](NATIVE_GS.en.md#scf-disk) |
| FDTD with LD/electronic states; continue in another job | [protected coupled restore](#coupled-restore) |
| Existing loop; save and exit before wall time expires | [CheckpointSchedule](#timed-save) |
| User-written code; bind its own storage format | Capability and execution interfaces below |

Postprocessing output archives do not replace restart checkpoints. No automatic checkpoint is written unless explicitly requested.

## Current mixed-session checkpoint formats

Direct persistence uses `salmonpy.native-direct/4`; replay persistence uses `salmonpy.native-replay/2`. Older schemas, automatic conversion and the unidentified-old-checkpoint option are removed. Custom hosts may use unidentified states in the current format, but identity presence must match between writer and reader.

[日本語](CHECKPOINTING.md)

`salmonpy.checkpointing` provides common input checking and execution for built-in and user
calculation implementations. `execute_checkpoint` calls explicitly bound SAVE/LOAD functions after input checking. It returns detached JSON
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
and runtime identity with `identity`. An empty identity is omitted from metadata; the saving and receiving declarations must agree. Checks compare the
explicit declarations; model identity and file integrity remain caller-owned.

## Execute calculation implementation-owned persistence

Bind zero-argument callbacks with `CheckpointTarget`, then call `execute_checkpoint`
with operation=save or restore. Paths and payload formats belong to the calculation implementation;
use closures or `functools.partial` to supply arguments. A complete file-based
example runs with `python -m salmonpy.examples.custom_checkpoint NEW_OUTPUT`.

The same entry can also transport an in-memory calculation implementation token:

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
and a `states` mapping keyed by string point IDs. Opaque tokens, arrays and validation records
remain local; they are not copied or sent through MPI. The calculation implementation controls their
ownership. Callback return values, including None, are accepted; callbacks must
raise exceptions to report failure. A return value of False is not an error signal.

Capability and callback availability are checked collectively before the first
callback. Groups execute their own point order in synchronized slots, so groups
with different point counts are supported. An ordinary callback exception becomes
RuntimeFailure on every rank; later slots stop. Completed state changes remain.
A collective calculation implementation must complete its own MPI protocol before raising; the common
entry cannot unblock ranks stranded inside a calculation implementation collective or recover a
Fortran STOP/process crash. Hosts decide whether to discard affected candidates
or abort their solver job. MPI initialization, splitting and finalization stay
with the calling program. Metadata/hash authentication and completion publication also remain
calling program/calculation implementation policies; the solver sample preserves its existing inventory/marker.

## Serial user calculation implementation

The packaged example contains a complete user model, JSON persistence, payload
SHA256 verification, independent states and comparison with continuous execution:

```sh
python -m salmonpy.examples.custom_checkpoint NEW_OUTPUT
```

It saves at17 and advances only15 more steps through32; current is exactly equal
to continuous execution. Changing the model is rejected before its loader runs.
The recurrence is synthetic interface data, not a physical material prediction.

This is the input checking part of that workflow:

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
The calling program must authenticate payloads and check that the separate state prepared for restoration is fresh
before binding its own loader to the execution entry. The sample above illustrates metadata input checking;
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
preservation of the live states. discarding the prepared restore state never authorizes changing live
objects. Under the calculation implementation isolation contract, no original is passed to discard.

`execute_checkpoint(..., operation='restore', atomic=True, live=live)` is the same
protected path. `execute_checkpoint` defaults to `atomic=False`, which does not preserve original states on failure.
Use `restore_checkpoint` for protected restoration. Direct-only targets are rejected before LOAD in that path. LD/SBE states can be staged via new_state plus
snapshot restore; custom JSON factories are shown in the runnable example.

### solver numerical candidates (experimental)

Explicitly rebuilt FDTD called directly from Python, standard TDDFT and HSE expose
checkpoint_protocol=1: modes0 save,1 immediate restore,2 stage a detached
numerical candidate,3 publish prepared memory,4 discard. One candidate per solver
lifetime is allowed. Pause evolution, finish mode2 and agree success across every
calling program group before mode3. On handled failure every group discards with mode4 and
can continue its original. HSE retains validated exchange/history data in memory;
mode3 never reopens the checkpoint. Reject older libraries via the protocol query.

The internal `salmonpy._native_checkpoint.NativeRestoreCandidate` helper and
`tools/test_native_restore_candidates_mpi.py` demonstrate the common entry:
bind stage/discard callbacks, use the common API for collective preparation,
then call publish in each group and install preprepared Python A/current/clocks.
The generic holder switch alone does not finish solver publication. This is an
internal experimental contract, not a stable public C ABI.

Only mutable numerical candidates are detached. Initialized MPI/Libxc/FFT resources
remain in the existing context; this does not enable concurrent solver lifetimes
or continued use of an old alias after publication. Publication uses prepared
memory but includes allocating assignments: allocation failure/process death is
outside the recoverable contract. The ordinary mpi_mixed_native integration is described below.
Exact-layout/input/hash input checking is still required.

This contract covers handled exceptions with healthy participating MPI ranks and
correct isolated factories. It does not recover process crashes, communicator loss,
Fortran STOP or a calculation implementation which mutates/shares its original state during staging.

## Existing mixed solver calculations

The solver mixed sample uses both common input checking and execution, retaining solver-specific
fresh-state and solver-ABI checks. Payload loading checks the format, input identity, MPI layout and file hashes.
Custom direct-state/1 sample adapters may expose `checkpoint_capability` as an
object/mapping. The declared checkpoint_format/checkpoint_modes/
checkpoint_decomposition attributes are also inspected. The sample's callback/clock/
A/current requirements are described in [COUPLED.en.md](COUPLED.en.md).

Current input checking requires the same topology and declared format/identity. Payload
conversion, changed-topology redistribution, generic archive I/O and atomic rollback
remain calculation implementation/calling program responsibilities. The execution entry calls supplied functions. The saving backend owns payload correctness; an immediate LOAD alone does not preserve original states on failure.

<a id="coupled-restore"></a>
## Protected restore in the ordinary mixed sample

`mpi_mixed_native --direct-resume PATH` now defaults to protected restore.
Initialized field/standard TDDFT/HSE require candidate-capable runtimes.
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
input checking errors are rejected before LOAD as well.

After success retrieve current states through session.states; old SBE handles are
retired. Release errors are separate coupled.checkpoint_cleanup_errors diagnostics,
not a failed restoration. solver resources still have one lifetime; allocation
failure during publication, process death and failed MPI remain outside recovery.

Older ABIs, process-isolated HSE workers and custom direct-only LOAD targets are
rejected before loading. The ordinary example always uses protected restoration. The common API still
allows explicit `execute_checkpoint(..., operation="restore", atomic=False)` for custom immediate-load operations;
that mode does not preserve original states on mid-load failure. There is no automatic fallback.

`tools/test_mixed_protected_restore_mpi.py` takes the ordinary sample arguments,
injects a rank failure after loading a later electronic candidate, then checks
original preservation (including two field-local LD materials), continued old
evolution and exact retry results for SBE plus internal-MPI TDDFT/HSE.
Source-tree tests can explicitly select a built LD/SBE library with `--kernel LIB`;
installed-wheel users do not need that override.

Live rewind does not truncate solver output files already written. Manage the
trajectory with the saved frames plus subsequent step returns, or resume in a
fresh MPI job/output directory to keep solver output streams separate.

Memory advice includes concurrent candidate storage; see [MEMORY.en.md](MEMORY.en.md). Exceeding the ordinary sample budget only warns and does not refuse protected restore.

## Mixed sample identity and restart

`python -m salmonpy.examples.mpi_mixed_native` preserves `--run-id`, logical point IDs, material IDs, material-definition and initial-state SHA256 digests, MPI ownership, and the LD cell-material mapping. The default `mixed-example` is a fixture label; provide an explicit run ID for each experiment. Output and temporary directories are excluded from definition fingerprints.

The same identity appears in `identity.json`, normal-run `result.json`, and the `point-response` archive. Identified direct checkpoints use `salmonpy.native-direct/4`; replay checkpoints use `salmonpy.native-replay/2`. Numerical payload formats and update order are unchanged. Direct restore checks identity and MPI ownership before loading numerical state; save checks agreement across all ranks.

Old arguments, checkpoint schemas and missing-field defaults are not supported. Use the current API and checkpoint writer. The mixed example always records material and initial-state identities.

## In-memory snapshots with point identity

`PointState.snapshot()` wraps the calculation implementation-owned snapshot with immutable `salmonpy.point-snapshot/1` identity. `state.restore(saved)` checks run, point, material-definition and initial-state identity plus package version before calling the calculation implementation. Cross-point/run and raw calculation implementation tokens are rejected. Even with explicit matching run/point IDs, the calculation implementation must still validate its model, time step and other settings.

This is process-local rewind, not branching, cloning, disk persistence or an MPI-wide transaction. Unsupported user-written calculation implementations raise `SnapshotError`. calculation implementation restore-failure guarantees remain the calculation implementation's responsibility; identity rejection invokes no calculation implementation restore callback. Raw-state snapshot APIs are unchanged.

`python -m salmonpy.examples.point_ld /new/output/directory` also verifies cross-point rejection and exact current continuation after same-point restore.

## TDDFT called directly from Python rejection before update (limited scope)

Newly prepared and built standard TDDFT called in the Python process exposes `sp_tddft_error_protocol() == 1`. Validation is collective on the attached material communicator:10 invalid/nonfinite input, stale A or missing input file;11 rank disagreement;12 unavailable lifecycle or repeated initialization;13 exhausted step capacity;14 pending separate state prepared for restoration. These codes return unanimously before numerical mutation. Initial input checking codes10/11 do not consume the lifetime; Python initialization can be retried with corrected input. Failed work directories remain for diagnosis.

Python maps only protocol-declared unanimous rejection to `InputError`. Old libraries, unknown codes or nonfinite results retain `RuntimeFailure`, with no inferred state guarantee. Every solver-group member must enter the same call. Invalid pointers, MPI/process failure, allocation errors, deep solver STOP and failures after numerical updates are outside this input checking-only contract; see the limited deep-recovery contract below. Actual initialization still consumes the one-process lifetime. Fast HSE and FDTD called directly from Python error contracts are unchanged.

This guarantee concerns an individual TDDFT state's pre-update rejection. A coupled session may already have advanced another component; it remains failed under the existing driver policy. This does not implement general coupled-step rollback.

## Direct rollback after deep TDDFT numerical errors

New standard TDDFT called directly from Python exposes `sp_tddft_recovery_protocol() == 1` and retains a per-step numerical payload in memory. Collective checks cover propagated wavefunctions(20), density(21), Hartree field(22), XC potential/energy(23), current/total energy(24), and a protected propagator error(25). Every material-group rank agrees before returning. Errors restore accepted orbitals, densities, fields, propagation history, XC state, counters and original overlap flags without propagation replay. Failed attempts add no ordinary RT output records.

A state returned by Python `NativeTDDFTBackend` raises `NumericalError` from `advance(A_old,A_new)` only for protocol-declared unanimous recovery. Its `A` and `steps` remain accepted; correct the cause and retry from that time. No automatic pulse weakening or dt changes are applied. Old libraries and unknown errors retain conservative `RuntimeFailure` handling without inferred recovery.

The following scope requires `execution_mode="checked"` and a supporting calculation library.

| Electronic state | Supported conditions |
|---|---|
| Unpolarized | Periodic fixed-ion middlepoint PZ, automatic/fixed TB-mBJ, Libxc-PZ/PBE |
| SOC / noncollinear spin | Periodic fixed-ion CPU PZ / Libxc-PZ/PBE, through a SOC-enabled interface such as retained-GS `gs.tddft_model(...)` |

If SOC SCF used an auxiliary spin field, set it to zero with `state.set_scf_spin_field([0., 0., 0.])` before starting TDDFT, then continue SCF without the field and inspect convergence, residuals and magnetization changes. `gs.tddft_model(...)` rejects a nonzero auxiliary field. Removing the field alone does not establish zero-field convergence. See [the auxiliary-field SCF procedure](NATIVE_GS.en.md#stabilize-soc-scf-with-an-auxiliary-spin-field).

Handled SOC numerical failures preserve wavefunctions, density, magnetization/spin-current and current output, allowing retry from the same completed step. CPU spatial SOC requires an orthogonal cell and `sp_tddft_soc_spatial_protocol()==1`. See [NICKEL_SOC.en.md](NICKEL_SOC.en.md) for tested conditions.
`execution_mode="performance"` does not promise this automatic checking or recovery. Use explicit guarded steps or restore a saved state when needed.
This guarantee does not cover +U, jellium, ionic motion/force output, Gram–Schmidt, projection/unfold or microscopic Maxwell coupling. Allow for backup/checking memory and runtime cost; budgets remain warning-only.

This is a TDDFT-group guarantee, not automatic rollback of an entire mixed FDTD step whose other components may have advanced. MPI/process death, allocation failures, uncontrolled third-party STOP, deep initialization STOP and output-device failures are excluded. One-rank fault injection exists only in a separately compiled `SALMONPY_RECOVERY_TESTING` build; production builds export no injection entry.

## Deep numerical recovery for HSE

New libraries for `NativeFastHSEBackend` declare `sp_fast_hse_recovery_protocol() == 1`. The guarded configuration is periodic, unpolarized, fixed-ion HSE06 with Taylor4 / learned PT-CN and internal k-only MPI1/2. Conventional HSE is not a fallback.

Each step backs up electronic arrays and ACE, factor/PT-CN histories, exchange caches and counters. FFT plans, MPI communicators and other owning resources are not copied. Collective checks cover finite orbitals, density, Hartree, XC and endpoints, plus handled ACE operator errors, PT-CN convergence/Gram/retraction checks and history updates. Handled failures directly restore accepted numerical state without replaying propagation or exchange history. A, clocks, current and energy remain unchanged; failed attempts append no ordinary RT records. Diagnostic logs may remain.

Only a declared protocol and unanimous known error permit Python to raise `NumericalError` and retry with corrected input at the same accepted time. Old libraries and unknown errors retain conservative `RuntimeFailure` behavior. Neither dt nor drive amplitude changes automatically. Per-step backups/checks add memory and runtime cost; budgets remain advisory.

The guarantee covers the HSE material group, not automatic rollback of an entire mixed step containing FDTD and other materials. Deep initialization, allocation failure, process/MPI termination, uncontrolled external STOP and output-device failure remain excluded. Do not infer recovery after an unhandled failure.

Run `tools/test_fast_hse_deep_recovery_mpi.py BUILD GS NEW_OUTPUT --steps 32 --fault-step 13` with two ranks, OMP/BLAS1 and an external timeout. `--split` uses two independent two-rank groups. `--overflow-A 1e100` stresses real PT-CN numerical divergence; it is not a scientific input. `--faults` requires a separate explicit `SALMONPY_RECOVERY_TESTING` build; production exports no injection entry point. See `docs/fast-hse-deep-recovery-validation-2026-11-17.md`. These small synthetic cases verify implementation, not converged physical accuracy.

## FDTD called directly from Python numerical recovery

New standard MPI managed by the calling program libraries declare `sp_eh_recovery_protocol() == 1`. `NativeEHField.advance` selects a checked entry and backs up field/PML, every field-local LD, deposited current and clock in memory. Collective gates check E after update (20), LD after update (21), H after update (22), and sampled fields/integrated A/ordinary observation arithmetic (23). Handled failures restore accepted state directly without propagation replay. A declared known recovery raises `NumericalError`, preserving Python A/E/clock/failed so the caller can correct input and retry. Failed attempts append no ordinary observation records; diagnostic/progress stdout may remain.

The guarded configuration is the existing FDTD updates of E and H plane-wave path with local LD. Dedicated linear response, ASE/ART and frequency-resolved observation are rejected before updating. Batch FDTD retains its old path. Old binaries without the checked entry and unknown failures remain conservative. Initialization, allocation, MPI/process termination, external STOP and output-device failure remain excluded. This protects the field group alone; it does not automatically roll back or retry a whole `CoupledStepper` / `ResponseSession` step. Budget for backup/check overhead; memory limits remain advisory.

Run `tools/test_eh_deep_recovery_mpi.py BUILD NEW_OUTPUT` with one/two ranks, OMP/BLAS1 and an external timeout. `--split` uses two independent two-rank groups. `--faults` requires a separate explicit `SALMONPY_RECOVERY_TESTING` build; production exports no injection entry. Overflow stress values are not scientific inputs.

See [NATIVE_STEPPING.en.md](NATIVE_STEPPING.en.md) for solver numerical/ordinary-output rollback after a step failure. Checkpoint LOAD is not a substitute for the step backup; concurrent use is rejected. Rebuild step targets when publication replaces owners.

## solver SCF memory snapshots through the common API

A live computed NativeGSState exposes `checkpoint_capability` and
`checkpoint_target(snapshot=None)`. These are collective on its internal group;
binding queries the runtime protocol and live SCF metadata, but performs no save or numerical iteration.
The format is `salmonpy.native-scf-memory/1`. Identity binds the material/model,
runtime and live calculation object token; even identical-model sibling owners differ.
Decomposition is the actual `(k, orbital, x, y, z)` tuple.

```python
from salmonpy.checkpointing import execute_checkpoint
comm = state.session.comm
plan = [{"name": "scf", "ranks": list(range(comm.size)), "point_ids": ["film"]}]
saved_run = execute_checkpoint(comm, plan, states={"film": state.checkpoint_target()})
with saved_run["states"]["film"] as payload:
    state.continue_scf(steps=5)
    execute_checkpoint(comm, plan, operation="restore", saved=saved_run["compatibility"],
                       states={"film": state.checkpoint_target(payload)})
    state.continue_scf(steps=5)
```

Only compatibility metadata is exchanged. NativeSCFSnapshot stays rank-local in
Fortran memory and must be explicitly closed. JSON metadata alone cannot restore
it, and it is not a process/disk restart. Ordinary PZ/Libxc-PZ/PBE and full-k CPU
HSE06 retain their existing scope and numerical implementation.

`execute_checkpoint(operation="restore")` remains **direct execution**.
Older runtimes without SCF candidate protocol1 provide no stage_restore, so
`restore_checkpoint` rejects those direct-only targets before numerical callbacks.
New matching runtimes support the staged path below. A common execute_checkpoint restore can leave earlier components
changed if a later callback fails; it does not promise group-wide rollback.
Wrong-owner/closed tokens and incompatible declarations reject before their
numerical callbacks. Do not substitute destructive restore for candidate staging.
`python -m salmonpy.examples.scf_resume BUILD_JSON NEW_OUTPUT` uses this entry;
add `--xc hse06` with the matching HSE build.

### Protected SCF staging and explicit solver publication

A matching runtime with `sp_gs_scf_candidate_protocol()==1` adds independent
**numerical** candidates to `state.checkpoint_target(payload)`. These are not
new runnable GS solvers: FFT/XC/MPI resources remain in the live original owner.
Staging deep-copies the snapshot into a detached Fortran payload, performs no
SCF initialization/orbital solve/exchange calculation, and reserves that owner.
Metadata can still be read; evolution/direct restore/snapshot/close are rejected
until the candidate is published or discarded. Calls remain collective and on
the main thread, with the original memory-only/profile/topology restrictions.

```python
from salmonpy.checkpointing import CheckpointState, restore_checkpoint
holder = CheckpointState(states={"film": state})
prepared = restore_checkpoint(
    comm, plan, live=holder, saved=saved_run["compatibility"],
    states={"film": state.checkpoint_target(payload)})
# All participants are now ready. Until this point solver originals are unchanged.
for candidate in prepared["states"].values():
    candidate.publish()  # collective within its own native group; returns state
state.continue_scf(steps=5)
```

The holder contains readiness candidates, not independently runnable states.
Publish every solver candidate in the same group order before resuming work.
A handled staging failure invokes discard callbacks (also for partly staged
slots), clears reservations and leaves original numerical data intact, so the
original trajectory can continue and restoration can be retried with new targets.
Before the ready vote, call candidate.discard() to abandon a manually staged
candidate; do not continue a partially published bundle after publication errors.
solver publication uses numerical assignments and may allocate. Fatal allocation,
process/communication failures, and arbitrary errors during solver publication
are outside this handled staging-error guarantee. This is not general rollback.
No files or numerical propagation are used by publish/discard. Snapshots must
remain open through staging/publication, and are closed explicitly afterwards.

Use `python -m salmonpy.examples.scf_resume BUILD_JSON NEW_OUTPUT --protected`;
add `--xc hse06` for a matching full-k HSE build. Old runtimes keep direct
mode; there is no destructive fallback when protected staging is unavailable.

solver SCF same-build/input/decomposition disk checkpoints bind through
`state.checkpoint_target(directory=...)`. scf_disk_capability is metadata-only;
staging authenticates SHA256 and inputs. Publish every returned candidate after
the ready vote, and recreate targets before retrying a failed restore.
See [NATIVE_GS.en.md](NATIVE_GS.en.md) for a two-state example and scope.

Prepared solver SCF candidates support close/context-manager discard. GS-owner
close discards its Python-owned candidate; snapshot close discards a candidate
derived from that snapshot. Direct C-calling program reservations are not automatically
discarded. Matching MPI call order is still required.

<a id="timed-save"></a>
## Periodic checkpoints and elapsed-time save-and-stop

Automatic checkpoints are disabled by default. Their cadence is separate from current, magnetization and other observation output. `CheckpointSchedule` calls an explicitly supplied save function at accepted step boundaries.

```python
from salmonpy.checkpointing import CheckpointSchedule

# Use the communicator of the entire coupled calculation; omit it for serial use.
# Construct before initialization if initialization must count against the budget.
schedule = CheckpointSchedule(
    interval_steps=1000,       # optional: save and continue every 1000 steps
    interval_seconds=600,      # optional: save and continue every 600 wall seconds
    stop_after_seconds=3300,   # optional: save and stop after 3300 wall seconds
    initial_step=0, comm=comm)

for step in range(restored_step + 1, final_step + 1):
    advance_one_step()
    if schedule.after_step(step, save_all_states):
        break
```

The caller supplies `advance_one_step` and `save_all_states(step)`. Use new destinations named by step. Coupled calculations must save the field, LD and electronic states together through their existing bundle-save interface. Saving only the electronic state is insufficient. Storage format remains backend-owned.

- All options default to None: no automatic save or extra loop MPI communication.
- When multiple intervals are specified, reaching any one triggers one save. Intervals count from the last successful save.
- Time means wall seconds since construction, not simulated physical time. Rank0's monotonic clock decides; every rank calls at the same completed step. An in-progress step is allowed to finish.
- A True return means saving succeeded and the caller must leave the loop and close resources. The library never terminates the process or finalizes caller-owned MPI.
- A save failure is reported to all ranks, not treated as successful completion. Retry at the same step; follow the backend's partial-file rules and use a new destination if needed.
- Reserve time before the scheduler deadline for one full step, checkpoint I/O and cleanup. Work before constructing the schedule is excluded.

The `mpi_mixed_native`, `mixed_hse_models` and `molecular_feedback` examples accept `--checkpoint-interval-steps N`, `--checkpoint-interval-seconds SECONDS` and `--stop-after-seconds SECONDS`. They save internal states directly into new directories such as `checkpoint-00001000`. Existing step-number save-and-stop options remain available.

Output joining was verified against actual fresh-process Ni SOC disk-restart files at1 and2 ranks: all8 current samples and the2 available magnetization/spin-current times exactly match continuous execution. This checks software restart wiring, not physical convergence.

### Save by elapsed time, then resume in another job

This example uses prepared calculation programs and a GS bundle. Four ranks are the default placement: one field rank, one SBE rank and two TDDFT ranks. Replace input paths with your prepared files and run on a compute node.

```sh
mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native \
  BUILD_JSON RUNTIME GS_BUNDLE first-job --steps 128 \
  --checkpoint-interval-steps 50 --stop-after-seconds 3300
```

Each destination includes the completed step, such as `first-job/checkpoint-00000050`. On a timed stop, the actual saved path is printed. Substitute that path for `SAVED_CHECKPOINT_DIRECTORY` below; do not predict the step at which the wall-time limit will be reached.

```sh
mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native \
  BUILD_JSON RUNTIME GS_BUNDLE second-job --steps 128 \
  --direct-resume SAVED_CHECKPOINT_DIRECTORY --stop-after-seconds 3300
```

Keep inputs, GS, calculation programs, total step capacity and MPI placement identical; choose a new output directory. This sample uses `--direct-resume` for timed direct checkpoints; its `--resume` is a separate history-replay route. The molecular and `mixed_hse_models` samples instead use `--resume` for their direct checkpoints. Check the chosen sample's arguments.

If all steps finish before the elapsed-time deadline, no timed final checkpoint is made. Periodic saving occurs only when a requested interval is reached. This is a coupling-verification sample; normal completion also performs independent reference calculations. Include that work in runtime estimates.

`after_step` returns True only for a successful save-and-stop. It returns False when no save is due or a periodic save succeeds and execution continues. Save failures raise an exception.
