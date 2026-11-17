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
