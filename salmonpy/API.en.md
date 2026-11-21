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

`fast_hse_mpi.NativeFastHSEBackend`: experimental host-MPI fast-PR HSE, k-only one/two-rank material groups and same-topology direct restart. Legacy create_state has one native lifetime/process. open_model returns NativeFastHSEModel for independent same-model full-k states; see SHARED_HSE.en.md and COUPLED.en.md.

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

For explicit sequential heterogeneous standard-model switching in one process, see [TDDFT_REGISTRY.en.md](TDDFT_REGISTRY.en.md).

[Retained GS handles and Python](NATIVE_GS.en.md).

`salmonpy.gs_input.GSInputMetadata` / `describe_gs_input` and retained
`state.input_metadata()` expose immutable declared-input descriptions. The optional
metadata extra is loaded only on explicit use; source fingerprints do not certify
physical GS identity or native restore compatibility. See the retained-GS guide.

`salmonpy.preparation.preflight_inputs(directory, manifest_sha256=None, registry=None)`
returns a JSON-ready input-admission report without starting a solver. Explicit
trusted resolvers are optional; without one the scope is file integrity only.
See PREPARATION.en.md for batch-job exit codes and per-rank startup reports.

## Declared native execution contracts

```python
from salmonpy.api import native_contract
contract = native_contract('fast_hse', route='shared')
assert contract['model_settings_owner'] == 'explicit model'
assert contract['runtime_verified'] is False
```

`native_contract` describes built-in `tddft`/`fast_hse` routes `legacy`, `shared`,
`registry`, and `fdtd` route `legacy`. It loads no library, executes no calculation
and enforces no new restriction. The result distinguishes process-compatibility
settings from explicit model ownership, caller-owned MPI, sequential main-thread
operations, entry-point collective ordering, current sign and failure scope.
Legacy standard-TDDFT C current is matter J=-j; its Python adapter returns charge
j. Fast-HSE and standard handle C currents are charge j. Follow exact headers for
query locality; a local diagnostic is not automatically collective.

This is declared documentation, not a live capability check or a stable-C-ABI
promise: `runtime_verified` is always false. Query actual runtime/checkpoint
capabilities separately. Unknown/custom component or route names return
`known=False`, `status='unverified'` and do not block a custom backend; that
backend owns its contract. Each result is independent and may be annotated by
the host without changing later reports. Existing calls, numerical state,
legacy fallback and accepted calculations are unchanged. Stable C/Fortran host
ABI conformance remains a separate milestone.

### Explicit live protocol probe

```python
from salmonpy.doctor import probe_native_protocols
report = probe_native_protocols('audited-build.json', component='tddft', route='shared')
# Inspect errors, observed protocols and protocol_match; runtime_ready stays unverified.
```

CLI: `python -m salmonpy.doctor --native-build audited-build.json --component tddft --route shared --json`.
Unlike `native_contract` and the separate static checker, this explicitly loads
native code in a child process after checking the recorded library SHA256 (and
rechecking in the child). A hash proves identity, not trust or build-audit
compliance. Use only trusted, audited libraries: loading can run constructors.
The probe calls only known zero-argument version functions; it creates no solver
state and does not request MPI initialization/attachment, compilation, installation
or calculation. Unknown/custom routes do not even open the record or load code.

`protocol_match` describes the selected version queries only. Missing/newer
versions are reported, never converted to automatic backend rejection or feature
fallback. `runtime_ready` remains `unverified`, and operation signatures, physical
profiles and checkpoint/topology compatibility are not certified. A child failure
is reported without taking down the host; `library_loaded=None` means the failed
child could not report whether loading completed. Keep this opt-in diagnostic
separate from `check_job`, which never loads the supplied code.

Query-only C and ISO_C_BINDING Fortran hosts are in
`tests/native/test_native_protocol.c` and the `SALMONPY_ABI_PROTOCOL_PROBE` mode of
`tests/native/test_io_ownership.f90`. Development runner:
`tools/test_io_ownership.py BUILD_RECORD NEW_OUTPUT --abi-contract` (add
`--fast-hse` for a fast-HSE library). It checks runtime SHA, audits every linked
application object, compares returned versions and verifies MPI stays uninitialized.
This is version-query interoperability, not full model/state/step ABI conformance;
the native ABI remains experimental.

### Explicit state-lifecycle host test

`tools/test_io_ownership.py BUILD_RECORD NEW_BUILD_DIR --abi-state --fixture EXISTING_PZ_FIXTURE --state-output NEW_RESULTS_DIR --ranks 2` runs C and ISO_C_BINDING Fortran hosts against an audited standard-TDDFT library. Add `--fast-hse` with a fast-PR-HSE build record and an existing full-k HSE fixture to exercise the same checks. The HSE fixture explicitly sets k-only MPI decomposition for each host run; it does not modify the supplied GS. Add `--state-symmetry` for the reduced-k fixture: it references the supplied `sym.dat`, enables the symmetry-aware HSE initializer and drives z, the direction preserved by this fixture's retained operations. The runner checks that representative k points are actually reported. This is a fixture-specific direction, not a restriction on full-k x/y/z driving. The caller initializes/finalizes MPI. The hosts compare two live states with independent runs, reject a mismatched previous A without changing accepted state, reject stale handles after slot reuse, and confirm a fresh state starts independently. Existing GS files are read in place and checked unchanged; this test neither copies GS data nor writes electronic checkpoints. Select results storage explicitly.

The common directory-scope declarations are in `native_context.h`, included by both handle headers. On macOS, use the existing `tools/build_embedded_tddft.py::mpi_link_command` policy when linking MPI hosts too: private executable copies of OpenMPI Fortran COMMON sentinels can corrupt MPI_IN_PLACE operations across the shared-library boundary. The runner applies this policy and records link commands and executable hashes. Do not infer correctness from a query-only probe.

Verified scope: standard PZ and fast-PR HSE full/reduced-k fixtures, C/Fortran hosts, MPI 1/2 ranks, paired versus independent currents. HSE uses its own documented signatures and rejection status (1 for mismatched previous A, versus 10 in standard TDDFT). The adapters live only in the common host test; they do not change either public ABI or numerical solver. The reduced fixture runs symmetry-aware PT-CN with learned exchange, 8 full k points represented by4 points. Other XC profiles, geometries and platforms require their own tests; the public native ABI remains experimental. These small fixtures check interoperability and state ownership, not converged material properties.

### Prepare the PZ host fixture

The standard `--abi-state` fixture needs `rt.inp`, `data_for_restart/` and `pseudo1_rps.dat`. Use `theory='tddft_pulse'`, periodic fixed ions, unpolarized spin, `middlepoint` and `trans_longi='tr'`. Omit `&parallel` from the RT fixture: the runner exercises both 1 and 2 ranks, so its generated template does not pin a decomposition. Cell, atoms, grid, k points, electron count, pseudo and XC must match the GS.

```sh
python tools/prepare_pz_state_fixture.py GS_INPUT EXISTING_RESTART PSEUDO NEW_FIXTURE
python tools/test_io_ownership.py BUILD_JSON NEW_BUILD --abi-state --fixture NEW_FIXTURE --state-output NEW_RESULTS --ranks 2
```

Use a Python environment with SALMONpy installed. `GS_INPUT` is the original GS input or a validated matching RT template. This periodic single-species PZ fixture preparer retains inline atoms and references existing restart/pseudo files: it never runs SCF, copies GS or downloads files. `examples/pz_host_gs.inp` is a small one-Si-atom starter with an 8³ grid and 2×1×1 k points. If needed, the user explicitly prepares GS on a compute node first. It is not a converged-material example.

`BUILD_JSON` may be a record file or a directory containing `build.json`. The GS-free `--parser-host` check verifies that status-only and explicit `native_host=.false.` calls retain ordinary checks; only explicit true selects field-specific checks.

`SalmonExecutable.from_runtime(directory)` delegates to `runtime.load_runtime`, verifying every bundle file and capability record. Prefer this entry over assembling a bare binary/hash manually.

## Optical conductivity and dielectric response: one polarization

Use `analysis.optical_conductivity(time, charge_current_x, total_E_x,
current_convention="charge", field_kind="total", units="atomic_gaussian",
window="hann", drive_rtol=1e-8, metadata={"material_id":"Si","polarization":"x"})`.
Then call `analysis.dielectric_function(sigma)`. Results have omega, values,
valid and config. Inputs are simultaneous scalar time series; A is not E.
For A-driven crystals the caller must provide E with consistent differentiation
and time placement. Convert SALMON matter current J explicitly using j=-J.
Only Gaussian atomic units are accepted. One selected drive gives a single
polarization response, not a complete tensor; one circularly polarized run does
not determine all tensor elements.

The existing compiled NumPy FFT uses exp(+i omega t), normalized by dt. With
exp(-i omega t), sigma=j/E and epsilon=1+4*pi*i*sigma/omega. Fortran computes
ratios and diagnostics. Drive bins at or below drive_rtol times the maximum,
and nonfinite conversion results, carry NaN with valid=False. Dielectric DC is
undefined. An all-zero drive yields masked output rather than stopping a solver.
External-field ratios are not silently treated as bulk dielectric functions.
The same window is applied to both signals; window, time range, dt, peak E and
caller provenance are retained. Windowed finite-time ratios are estimates;
check duration/window dependence, especially for broadband drives.

`optical_linearity(reference, candidate, rtol=...)` compares supplied normalized
conductivities on common valid frequencies, reporting relative RMS. The caller
must ensure matching material, polarization and drive shape. Zero reference
response has undefined relative error/acceptance. `optical_sum_rule(sigma,
expected_weight=..., normalization="...")` integrates Re sigma between adjacent
valid bins with the trapezoidal rule. For an all-electron Gaussian-AU system,
the positive-frequency f-sum reference is pi*n/2. Pseudopotentials, nonlocal
current operators and missing bands require an appropriate caller-specified
reference. Masked gaps and unobserved frequency ranges are not extrapolated.
Neither diagnostic is a physical accuracy certificate or launches another run.

Run `python -m salmonpy.examples.optical_response` for synthetic data, not a
converged real-material optical calculation.


Additive API: salmonpy.preparation.prepare_fdtd(calculation, directory, mpi=None) writes a new native FDTD input folder and returns its Path only. FDTDMPI.rank_arguments and command(executable) configure/preview argv, without execution. See PREPARATION.en.md.
