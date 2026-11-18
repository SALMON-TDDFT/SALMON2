# Public Python API and compatibility

[日本語](API.md) · Topic: `public_api` · Paired revision: `2026-11-14.3`

Version **0.1.0a5** exposes `salmonpy.__version__`. The distribution metadata has
the same version. This is a preview, not a stable 1.0 API. The material C ABI preview1 is documented in [C_ABI.en.md](C_ABI.en.md).

## Public boundary

`salmonpy.api.PUBLIC_API` is the explicit module-to-name inventory. Import from
the named module; incidental imports and underscore-prefixed implementation
names are not public. For returned states/results, the methods and properties
specified in their topic guides are public; direct construction of private
field states is not. `StateSnapshot` is a token type, not a caller constructor.

| Modules | Intended entry |
|---|---|
| `ld`, `sbe`, `crystal`, `materials` | Materials, initial states and A-driven response |
| `maxwell`, `timing` | Field factories, coupling and time declarations |
| `fields`, `geometry`, `analysis`, `outputs` | Pulses, structures, analysis and archives |
| `dft`, `bands`, `tddft`, `fdtd`, `runtime` | Explicit external solver execution and setup |
| `ase`, `ase_calculator`, `salmon_vr` | Optional integrations |
| `distributed`, `material_response`, `embedded_tddft`, `coupled` | Experimental MPI/embedded APIs with their existing lifetime restrictions |
| `errors`, `snapshots`, `api` | Error families, opaque state tokens and inventory |

The root `__all__` intentionally exports only the version, `LorentzDrude` and
error types. It does not eagerly import optional integrations. Import performs
no compilation, download, native loading or MPI initialization. No console output
is added by snapshot/restore. Numerical dependencies remain NumPy plus explicitly
prepared LD/SBE kernels; optional/native solver requirements have not changed.

## Compatibility policy

Every distributed behavior/API change gets a new package version. Within a
release, do not change units, current sign, time layout or array names silently.
Preview releases may introduce documented breaking changes: pin the tested
package version and consult the migration note. A future stable release should
use a deprecation period before removals; no such stability is promised for all
experimental MPI/embedded functions today. This inventory is an API boundary,
not a claim that every listed function has equal maturity.

Archive, timing, runtime and snapshot schemas have independent versions.
Do not infer archive compatibility from the package version alone. In-memory
snapshots require the exact package version and originating model object.

Migration to a5: numerical update formulas and current conventions are unchanged.
Crystal interfaces use **A/c → charge-current density j**, with native electronic
matter current J = -j. The earlier a4 development history changed this convention;
old data must be interpreted using its recorded metadata, never blindly negated.
Bulk TDDFT remains transverse by default. Longitudinal feedback is explicit.

## Error contract

| Type (`salmonpy.errors`) | Meaning | Legacy catch retained |
|---|---|---|
| `SalmonpyError` | Base for explicitly raised package library errors | `Exception` |
| `InputError` | Invalid values, compatibility or state transitions | `ValueError` |
| `InputTypeError` | Unsupported input type; also an `InputError` | `TypeError` |
| `SnapshotError` | Invalid/incompatible restore; state unchanged | `ValueError` |
| `NumericalError` | Detected nonfinite/invalid numerical output | `FloatingPointError` |
| `RuntimeFailure` | Missing/failed runtime or result not established | `RuntimeError` |

Existing `RuntimeUnavailable` and `OutputSaveError` remain available and now
inherit `RuntimeFailure`. Catch standard families as before, or use the new
package families. OS/file errors (including `FileExistsError`), process timeouts,
NumPy conversions, missing optional dependencies and user callback errors may
retain their original Python types. We do not blanket-wrap external exceptions.
A missing convergence report is a runtime failure, not proof of numerical divergence.

An exception family does not imply rollback for a whole calculation. Built-in
LD/SBE/Python-field step and restore failures preserve the accepted state;
MPI response/coupling failures may invalidate the session. Native STOP/abort is
not converted into a recoverable Python error by this change.

Process-isolated `TDDFTState` also exposes `save_checkpoint(directory)` and
`restore_checkpoint(directory, max_bytes=None)`. Supply an identical initial
model explicitly to resume from Zarr history in another process. GS payloads are
not included. See [state persistence](SNAPSHOTS.en.md).

For one serial model contract, the explicit conformance checker and the point-MPI adapter, see [MATERIALS.en.md](MATERIALS.en.md).

`fast_hse.FastHSEMaterialProvider` / `FastHSEResponseBackend` are preview entries for an explicitly built serial fast-PR worker in outer-MPI response groups; see COUPLED.en.md. They do not initialize inner HSE MPI.

`fast_hse_mpi.NativeFastHSEBackend`: experimental host-MPI fast-PR HSE, k-only one/two-rank material groups and same-topology direct restart. One native lifetime per process, main thread only; see COUPLED.en.md.

`checkpointing.CheckpointCapability` and `check_checkpoint_compatibility` provide shared declaration/admission for user backends; see [CHECKPOINTING.en.md](CHECKPOINTING.en.md).

`checkpointing.CheckpointTarget` and `execute_checkpoint` invoke backend-owned persistence after admission and return opaque results locally.

`CheckpointState` / `restore_checkpoint` stage independent candidates before collective publication. AtomicRestoreError distinguishes preservation of live holders from candidate cleanup failures.

Advisory memory APIs are MemoryEstimate/MemoryReport/memory_advice/process_memory in salmonpy.memory. Excess only warns; this is not admission or an RSS cap. See [MEMORY.en.md](MEMORY.en.md).

`stepping.StepTarget`, `guarded_step` and `StepRecoveryError` are preview APIs for explicit complete backend rollback contracts. See the paired STEPPING guides.

`native_steps.native_step_targets` is the preview bridge from audited native state/ordinary-output transactions to common protected steps. See the paired NATIVE_STEPPING guides.

### Ion inputs during solver initialization

The standard solver and fast-PR HSE now capture atom/element counts, positions
(Cartesian or reduced) and species IDs into the internal value object
`s_ion_initial_config`. `initialize_ions` resolves against the selected system
lattice and creates independent owned arrays. Invalid dimensions, nonfinite
positions, out-of-range species and non-one-based input arrays reject before
changing the target. The explicit helper never reads/writes shared input.
The compatibility initialization path retains its old reduced-to-Cartesian input
normalization for downstream consumers, with unchanged valid calculation results.

This is an internal ownership boundary, not a new Python API or stable C ABI.
Cell/k/symmetry/communicator and other initialization dependencies remain; it does
not establish arbitrary different-input GS/TDDFT contexts in one process. Memory
allocation failures and global initialization STOP paths remain outside recovery.

### Explicit cell and real-space grid resolution

The standard solver now shares fast-PR HSE's value interfaces: `s_cell_config`,
`s_cell_geometry`, `s_rgrid_config`, `resolve_cell` and `init_grid_whole_config`.
Both profiles use capture-only `cell_config_from_input`; the explicit resolvers
read no global input. They preserve valid legacy cell precedence, arithmetic,
periodic/isolated index placement, epsilon rounding and Maxwell consistency
lengths. Rejected resolutions retain the target geometry/index arrays. These
helpers resolve cell values and whole-grid indices; they do not initialize all
coordinates, parallel decomposition, symmetry or FFT resources.

Existing initialization wrappers still capture legacy inputs and normalize `al`
for downstream consumers. This is not complete different-cell GS/TDDFT reentrancy
or a new public Python/C interface. The optional RT-step structure profile carries
the same value types; preparation verifies its new before/after source hashes.

### Explicit initial orbital layout

Standard and fast-PR initialization now capture mode, periodicity, spin/SOC,
state/electron counts, temperature and the k-mesh selection into internal
`s_orbital_initial_config`. `resolve_orbital_layout` chooses spin count, orbital
count and real/complex storage without reading shared inputs. Resolved cell
orthogonality is supplied explicitly; fast-PR uses its owned occupation settings,
effective k mesh and the same captured legacy hybrid predicate. Standard
file-supplied k points retain the legacy `-1` mesh sentinel for this decision.

Valid GS/RT selection rules, including the `temperature < -1e-12` RT branch, are
unchanged. Invalid mode, nonfinite temperature or an unusable resulting layout
reject without changing the result object. This is selection of initial array
layout, not full occupation-capacity/physical validation, SCF-input ownership or
recoverable native initialization. No new public Python/C API is added; remaining
k-point/symmetry/communicator and other solver dependencies still apply.

Fast-PR is the only supported HSE engine. `NativeFastHSEBackend(build_record, work_root,
allow_symmetry=True)` explicitly permits scalar periodic CPU/full-q/internal-k1/2,
Taylor/PT-CN symmetry with matching mode capability. Shared symmetry settings are immutable, one native
lifetime per process. The default False and existing C init remain symmetry-denied.
Internal-k2 requires sp_fast_hse_symmetry_kmpi_protocol()==1. Taylor strict/learned and PT-CN require sp_fast_hse_symmetry_history_protocol()==1; learned PT-CN keeps interval=1 and predictor_corrector=n.

`sp_fast_hse_init_symmetry` has the existing init arguments and declares
`sp_fast_hse_symmetry_protocol()==1`. Only unanimous finite status34 with that
protocol becomes a state-preserving Python `InputError`; correct A before retrying.
Unknown/old ABI errors do not imply recovery. Use `step_protected` for whole FDTD
loop retry. Common direct/staged checkpoint integration is available.

See [the symmetry guide](FAST_HSE_SYMMETRY.en.md) for mixed samples, exact scope and
stable near-zero-gradient XC evaluation and regression checks. Connectivity/restart tests are not
physical-convergence certificates.

Reduced-k symmetry now supports two internal k ranks and same-partition direct resume with the separate native k-MPI capability. Checkpoint repartitioning remains unsupported; see the symmetry guide.

NativeFastHSEBackend collectively admits periodic cell/grid/atom inputs before native lifetime. Repair source inputs and construct a new provider after InputError. Restart compatibility and general native initialization safety remain separate; see the symmetry guide.

Missing the SCF threshold emits RuntimeWarning and permits usable finite results
in downstream workflows. `GroundStateResult.converged` and the saved JSON
`converged` field expose convergence. Manifest `status='complete'` denotes
completed execution, not convergence.
