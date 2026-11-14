# TDDFT: different materials and dt in one process

## Scope and status

The user selected **design for different materials and dt simultaneously** on 2026-11-15. This supersedes the same-configuration restriction in `2026-11-15-tddft-serial-session.md`. The previous scope produced an unvalidated local prototype only; no library was published and no public API was added. Keep its sources/failing test/build logs as evidence outside the checkout. This document is a design, not an implementation or a numerical certification.

The target is multiple material models and multiple time-evolution states in one native process. A model describes material/initial electronic data; each state independently fixes dt, initial A and capacity. Initially calls are serial on the main thread. Different dt does not imply concurrent native execution. The host schedules calls at each state's declared time; this change does not add temporal interpolation or a new field-coupling integrator.

Keep the existing process-isolated TDDFT API, replay checkpoints and current conventions. Continue using an external pinned source plus opt-in patch/adapter, not a full SALMON source import. Existing solver binaries and production jobs are untouched.

## Alternatives considered

1. **Explicit configuration and ownership (selected).** Material configuration and all mutable state/workspaces are explicit arguments/owners. Different grid, k mesh, XC, pseudopotential and dt cannot silently affect a sibling. This requires staged changes through the TDDFT dependency graph but matches the existing material ABI and long-term component design.
2. Save/restore every global module before/after each call. This might shorten an initial demonstration, but missing even one allocatable, hidden cache, file unit or opaque handle causes cross-talk. Copying unowned external resources is unsafe. Do not use whole-module swapping as the final design.
3. Keep separate worker processes. Already supported and retained as the reference and robust fallback. This does not establish multiple contexts in a single native process. Repeatedly loading renamed shared libraries is also not a portable ownership solution.

## FDTD clock policy (confirmed 2026-11-15)

Different materials may use different fixed dt only when they are independent states or belong to different coupling loops. Within one current FDTD/TDDFT loop, **every material state uses the same dt as the field**, even when materials, grids or XC differ. Preserve the existing admission checks; reject dt mismatch before any coupled step. Integer ratios such as dt_TDDFT=dt_FDTD/2 do not imply supported subcycling.

For native Yee E/H, E is at n*dt_FDTD, A/current are associated with the field's half-step clock, and the field consumes the previous current with lag1. A phase offset is not a remedy for unequal dt. A future multirate adapter needs its own clock/ratio declaration, A sampling/interpolation, returned-current timing/deposition and conservation/convergence tests. It is excluded from this context-extraction implementation; the backend must not silently repeat/skip steps, resample A, average j or override a state's dt.

Acceptance includes two different material timing declarations sharing field dt, mismatch at either finer or coarser integer/noninteger ratios, and equal-dt but wrong-half-step rejection. Clock tests are interface tests; they do not certify a multirate physical scheme.

## Ownership and interfaces

| Object | Owns | Lifetime |
|---|---|---|
| `rt_material` | validated immutable material config, initial orbitals/occupations, grid/cell/k definitions, pseudo definitions and hashes, XC definition | referenced by states; Python model close prevents new states, existing states retain a lease |
| `rt_state` | material lease, fixed dt/time0/capacity, A/current/iteration parity, evolving orbitals, density, Vh history, predictor/corrector data, mutable nonlocal projectors, XC instances/workspace, Poisson workspace, communication requests/buffers, output handles, timers, error state | destroyed independently; never frees a sibling's resources |
| `rt_runtime` | handle registry, serial-call guard, backend capabilities, allocator/accounting and explicitly borrowed host services | process/library lifetime; does not own an application's MPI lifetime |
| parsing context | namelist units and buffers, coordinate/pseudo input staging, unit conversion records, validation diagnostics | local to material creation; produces an immutable config, not active global settings |

Separate `rt_material_config` from `rt_step_config`: grid/k/spin/XC/pseudo and physical conventions belong to the model; dt, capacity, observation policy and initial endpoint belong to the state. Physical constants remain read-only process constants. Validate model/state inputs once and pass configuration explicitly to kernels. Do not change a native global dt when switching states.

Initially retain independent state arrays/workspaces. Share only proven immutable material data with reference ownership. Avoid intrinsic derived-type assignment for XC handles or MPI resources. Copy initial orbitals deliberately; allocate a new workspace rather than copying a live communicator/request. The fast path need not reparse files or clone all model data every step.

### Existing material ABI

Use `salmonpy_material.h` ABI1 semantics when the native implementation satisfies them:

- `model_create(config_bytes, length, model_handle, error)` loads an explicitly supplied validated model description/GS reference. Config is versioned; it never contains host pointers or an implicit active model.
- `state_create(model_handle, dt, time0, initial_A, state_handle, error)` fixes state-specific integration settings. Fresh GS time0 initially remains zero, as in the existing TDDFT API.
- `state_step(state_handle, Anew, time, j, error)` changes only that state; `state_observe` is read-only.
- Handles are positive 64-bit opaque IDs; zero and stale handles fail. Destroying then recreating a state cannot recycle a still-observable stale ID. Check integer capacities before allocating/advancing counters.
- C `model_destroy` reports BUSY while states retain it. Python `model.close()` may defer native destroy through a lease, preserving existing loader semantics. State close is idempotent in Python; raw C stale-handle destruction returns HANDLE.
- Keep atomic units, Ac=A/c input and **charge current j=-J_native**. Convert the sign exactly once. Each state's `salmonpy.time/1` includes its own dt/integrator/lag. Different dt is not silently coerced to a common clock.

Do not label an experimental adapter ABI1-compliant until it satisfies the ABI's no-process-termination and failed-step atomicity requirements. A temporary diagnostic native interface may remain explicitly experimental while these obligations are being removed from transitive code.

## Concrete source extraction boundaries

The following findings come from the safe prepared runtime examined on 2026-11-15; lexical locations are evidence, not a full transitive reentrancy proof.

| Current owner/path | Required change |
|---|---|
| `io/salmon_global.f90` and `io/inputoutput.f90` | Split immutable material config, step config and parser-owned mutable arrays/units. `read_input` must fill a new config without reallocating another model's atom arrays. Keep CLI namelist loading as a frontend to this same config builder. |
| `rt/embedded_tddft.f90` singleton state/used/accepted_A | Replace by registry entries and explicit material/state owners; move endpoint, iteration and capacity to state. Remove one-lifetime assumptions only when their underlying resources are owned. |
| `rt/time_evolution_step.f90` imports `salmon_global` | Pass config with dt, propagator and supported observation flags; resolve transitive dependencies too. Preserve odd/even orbital buffers and Vh history exactly. |
| `rt/taylor.f90` imports `inputoutput:n_hamil` | Pass Taylor order/config explicitly rather than sharing an active input module. |
| `xc/salmon_xc.f90` module scratch arrays | Make scratch shape/spin/context-owned. Same dimensions are insufficient proof of independence. |
| `structures:s_xc_functional`, Libxc instances | Initialize and finalize XC handles for their explicit owner; never shallow-copy pointer wrappers. Support differing functionals only after independent lifetime/shape tests. |
| mutable k/A-dependent projector data in `s_pp_grid` | Per-state data; changing one A cannot update a sibling's projectors. Immutable pseudo source data may be model-shared. |
| `parallel/parallelization.f90`, attached communicator in `communication.f90` | Explicit parallel context; serial first. Later duplicate supplied MPI communicators per owner, free only owned children, never finalize host MPI. |
| `common/sendrecv_grid.f90` | Keep buffers and requests under each state's srg; no copied live persistent requests. |
| `misc/timer.f90` | Per-context timing workspace; process-wide immutable timer IDs are fine. A sibling's init/close cannot reset/stop another state's timers. |
| fixed output units (e.g.030), ofl units, Python cwd changes | Context-owned `newunit`/paths or explicit output sink. No process cwd change in the final native API. Disable unsupported observations by validation, not by silently discarding requested output. |
| transitive STOP/error stop and nonfinite paths | Return structured status with owned diagnostic; unwind only resources from the failed operation. Audit supported branches before exposing ABI1. |

Move fields in compilable groups rather than pretending that changing the top-level derived type removes all global dependencies. Track remaining accesses to state/config module data in the supported dependency closure. Preserve legacy CLI output/numerics via a frontend-owned context; do not rewrite historical snapshots.

## Errors and transactional state

Validate handles, A shape/finiteness (scalar IEEE predicates), capacity and caller settings before mutation. Rejection must preserve current, time and owned resources. A failed creation returns handle0 and releases only its allocations/units/XC objects.

For late numerical/allocation errors, compute a candidate step with an explicit commit boundary. Initially correctness takes priority: keep accepted buffers and a separately owned candidate/workspace; swap only after all required finite/collective checks pass. Include Vh/predictor histories, endpoints and other mutable fields in the transaction. Measure the extra copy/memory cost and then reduce copying through proven buffer ownership. Do not claim atomicity from restoring just A/current/iteration. Native STOP is not an error return; until removed on all supported paths the adapter remains experimental and cannot replace the process-isolated public backend.

A host call guard returns BUSY on recursive/concurrent calls in the initial serial implementation. Never serialize only at Python level while claiming that raw C is safe. No concurrent threads/GPU claims in this phase.

## Memory and persistence

Report model, accepted-state, candidate/workspace, communication, output and temporary I/O bytes separately. Derive allocation capacities from actual grid/k/band/spin/dtype and checked 64-bit products. Optional caller budgets are explicit; no machine-independent memory ceiling or silent radius/accuracy reduction. Allocator retention and peak RSS are measured separately from live owned bytes. Closing states must lower live allocation counts even if the system allocator retains pages.

The existing GS/2 archive remains the immutable model interchange format. Material identity includes all physical configuration and referenced hashes. dt belongs to each state and is part of its snapshot/checkpoint contract. Existing process replay checkpoints are not native wavefunction snapshots and must not silently load into a new backend. Native checkpoints require a separate schema/capability once the complete mutable state is explicit.

## Validation gates

Use synthetic fixtures only until real GS collection is authorized. They test ownership/interchange, not converged-material physical accuracy.

1. One-context baseline versus unchanged audited process worker: 2/3/6/8 k points and nonuniform weights, at least26 steps covering both orbital parities. Current rtol1e-12, atol1e-14; report exact differences instead of assuming bit equality across different builds.
2. Same material with dt=.001 and .002, different A histories and unequal step counts; interleave irregularly. Each trajectory must match its own independent worker; verify each reported time.
3. Different material configs: grid8³ versus12³, distinct cells/atomic coordinates/pseudo hashes where fixtures permit, different k/band counts; compare every state to its own reference. Include A→B→A creation and both destruction orders.
4. Different XC PZ/Libxc-PBE after XC ownership extraction; test independent XC finalization, grid/spin scratch shapes, partial creation failures. Extend to TB-mBJ only with its additional mutable payload covered.
5. Destroy/recreate either state/model, loader close with retained states, stale IDs, zero/one/exhausted capacity, NaN/Infinity in every input component, invalid config, failed allocations/I/O and late numerical failure. Verify surviving states and transactional rejection.
6. Memory: one/two/many states, repeated create/destroy, returned live bytes, peak RSS, handles/units/descriptors. Test finite values and scalar NaN/Infinity rejection. No fixed universal cap.
7. Source and installed-wheel suites, fresh source/object audits for every new executable/library and GNU/Linux toolchains. Fujitsu/MPI/GPU remain explicitly unverified until separately tested. No remote production changes.

## Delivery stages

A. Extract config/parser ownership and make the single-context reference path accept explicit config without numerical changes.
B. Make dt/config use explicit through the supported RT call graph; prove two dt values for the same material.
C. Extract XC, projector, Poisson, I/O/timer lifetimes and prove two different material/grid/XC models.
D. Remove reachable native-fatal paths and complete transactional errors; only then expose the backend through material ABI1/Python preview.
E. Validate reference comparisons, memory, distribution and documentation; keep process isolation available. MPI ownership and native checkpoint design are subsequent milestones, not inferred from serial success.
