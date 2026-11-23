# Public Python API and compatibility

## Common query before C/Fortran calls

`sp_host_abi(capacity, values)` reports the current host contract and primitive
widths. Include `native_context.h` and supply at least10 `int64_t` words. Success0
writes only the first10; short capacity returns13 without writing. The query
allocates no workspace, communicates nothing and initializes neither MPI nor a solver.

| C index | Meaning | Current value |
|---|---|---|
|0|Host contract version|1|
|1–5|C int, C double, int64, Fortran INTEGER, real(kind=8) bits|32,64,64,32,64|
|6|Communicator encoding|1: MPI_Comm_c2f result widened to int64; never raw MPI_Comm|
|7|MPI initialization/finalization owner|1: caller|
|8|State-operation threading|1: main thread, sequential calls|
|9|Directory protocol required for I/O|3; not proof that I/O operations exist|

Fortran indexing starts at1. Matching this descriptor does not certify operation
signatures, MPI implementation, physical profiles or checkpoint compatibility;
query those separately. Built-in Python MPI attachment validates it after SHA256
checking and before attachment. A mismatch on one rank rejects the entire owning
group before attachment. There is no old-version fallback or implicit build.
User-written backends retain their own contracts instead of being forced into this one.

`salmonpy.api.query_host_abi(library)` queries an explicitly loaded trusted library;
mismatch raises `RuntimeFailure`. Built-in collective MPI admission reports an
`InputError`. To query in a separate process without running a calculation:

```sh
python -m salmonpy.doctor --native-build BUILD_JSON --component fdtd --route single
```

`single` means one state, `shared` means multiple states of one model, and
`registry` means registered different-input models. The old `legacy` route is not an alias.


[日本語](API.md)

This page describes the available functions and classes, errors and saved data.
To choose a calculation, start with [getting started](README.en.md) or the
[examples](examples/README.en.md).

Read the installed version from `salmonpy.__version__`. APIs can change during
this unreleased development period; use the matching package and documentation.
For calls from C or Fortran, see [C_ABI.en.md](C_ABI.en.md).

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
| `distributed`, `material_response`, `embedded_tddft`, `coupled` | Experimental APIs for calculations in the Python process and MPI work sharing |
| `errors`, `snapshots`, `api` | Error families, objects representing saved states and inventory |

The root `__all__` intentionally exports only the version, `LorentzDrude` and
error types. It does not eagerly import optional integrations. Import performs
no compilation, download, solver loading or MPI initialization. No console output
is added by snapshot/restore. Numerical dependencies remain NumPy plus explicitly
prepared LD/SBE calculation routines; requirements for optional packages and calculation programs have not changed.

## Current API policy

SALMONpy is unpublished and does not provide backward-compatibility branches,
old argument aliases, format migration or a deprecation period. Use the current
documented API, runtime and writers together. Units, current sign and time layout
must remain explicit. Package versions identify releases; pin the tested version.

Archive, timing, runtime and snapshot schemas have independent versions.
Do not infer archive compatibility from the package version alone. In-memory
snapshots require the exact package version and originating model object.

Migration to a5: numerical update formulas and current conventions are unchanged.
Crystal interfaces use **A/c → charge-current density j**, with solver electronic
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
MPI response/coupling failures may invalidate the session. solver STOP/abort is
not converted into a recoverable Python error by this change.

Process-isolated `TDDFTState` also exposes `save_checkpoint(directory)` and
`restore_checkpoint(directory, max_bytes=None)`. Supply an identical initial
model explicitly to resume from Zarr history in another process. GS payloads are
not included. See [state persistence](SNAPSHOTS.en.md).

For one serial model contract, the explicit conformance checker and the point-MPI adapter, see [MATERIALS.en.md](MATERIALS.en.md).

`fast_hse.FastHSEMaterialProvider` / `FastHSEResponseBackend` are preview entries for an explicitly built serial HSE worker in outer-MPI response groups; see COUPLED.en.md. They do not initialize inner HSE MPI.

`fast_hse_mpi.NativeFastHSEBackend`: experimental MPI managed by the calling program HSE, k-only one/two-rank material groups and same-topology direct restart. Legacy create_state has one solver initialization per process. open_model returns NativeFastHSEModel for independent same-model full-k states; see SHARED_HSE.en.md and COUPLED.en.md.

`checkpointing.CheckpointCapability` and `check_checkpoint_compatibility` provide shared declaration and compatibility checking for user calculation implementations; see [CHECKPOINTING.en.md](CHECKPOINTING.en.md).

`checkpointing.CheckpointTarget` and `execute_checkpoint` invoke calculation implementation-owned persistence after input checking and return opaque results locally.

`CheckpointState` / `restore_checkpoint` stage independent candidates before collective publication. AtomicRestoreError distinguishes preservation of live holders from discarding the prepared restore state failures.

Advisory memory APIs are MemoryEstimate/MemoryReport/memory_advice/process_memory in salmonpy.memory. Excess only warns; this is not input checking or an RSS cap. See [MEMORY.en.md](MEMORY.en.md).

`stepping.StepTarget`, `guarded_step` and `StepRecoveryError` are preview APIs for explicit complete calculation implementation rollback contracts. See the paired STEPPING guides.

`native_steps.native_step_targets` is the preview bridge from audited solver state/ordinary-output transactions to common protected steps. See the paired NATIVE_STEPPING guides.

### Ion inputs during solver initialization

The standard solver and HSE now capture atom/element counts, positions
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

The standard solver now shares HSE's value interfaces: `s_cell_config`,
`s_cell_geometry`, `s_rgrid_config`, `resolve_cell` and `init_grid_whole_config`.
Both configurations use capture-only `cell_config_from_input`; the explicit resolvers
read no global input. They preserve valid legacy cell precedence, arithmetic,
periodic/isolated index placement, epsilon rounding and Maxwell consistency
lengths. Rejected resolutions retain the target geometry/index arrays. These
helpers resolve cell values and whole-grid indices; they do not initialize all
coordinates, parallel decomposition, symmetry or FFT resources.

Existing initialization wrappers still capture legacy inputs and normalize `al`
for downstream consumers. This is not complete different-cell GS/TDDFT reentrancy
or a new public Python/C interface. The optional RT-step structure configuration carries
the same value types; preparation verifies its new before/after source hashes.

### Explicit initial orbital layout

Standard and HSE initialization now capture mode, periodicity, spin/SOC,
state/electron counts, temperature and the k-mesh selection into internal
`s_orbital_initial_config`. `resolve_orbital_layout` chooses spin count, orbital
count and real/complex storage without reading shared inputs. Resolved cell
orthogonality is supplied explicitly; HSE uses its owned occupation settings,
effective k mesh and the same captured legacy hybrid predicate. Standard
file-supplied k points retain the legacy `-1` mesh sentinel for this decision.

Valid GS/RT selection rules, including the `temperature < -1e-12` RT branch, are
unchanged. Invalid mode, nonfinite temperature or an unusable resulting layout
reject without changing the result object. This is selection of initial array
layout, not full occupation-capacity/physical validation, SCF-input ownership or
recoverable solver initialization. No new public Python/C API is added; remaining
k-point/symmetry/communicator and other solver dependencies still apply.

HSE is the only supported HSE engine. `NativeFastHSEBackend(build_record, work_root,
allow_symmetry=True)` explicitly permits scalar periodic CPU/full-q/internal-k1/2,
Taylor/PT-CN symmetry with matching mode capability. Shared symmetry settings are immutable, one solver
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

Reduced-k symmetry now supports two internal k ranks and same-partition direct resume with the separate solver k-MPI capability. Checkpoint repartitioning remains unsupported; see the symmetry guide.

NativeFastHSEBackend collectively admits periodic cell/grid/atom inputs before solver initialization. Repair source inputs and construct a new provider after InputError. Restart compatibility and general solver initialization safety remain separate; see the symmetry guide.

Missing the SCF threshold emits RuntimeWarning and permits usable finite results
in downstream workflows. `GroundStateResult.converged` and the saved JSON
`converged` field expose convergence. Manifest `status='complete'` denotes
completed execution, not convergence.

For explicit sequential heterogeneous standard-model switching in one process, see [TDDFT_REGISTRY.en.md](TDDFT_REGISTRY.en.md).

[Retained GS handles and Python](NATIVE_GS.en.md).

`salmonpy.gs_input.GSInputMetadata` / `describe_gs_input` and retained
`state.input_metadata()` expose immutable declared-input descriptions. The optional
metadata extra is loaded only on explicit use; source fingerprints do not certify
physical GS identity or solver restore compatibility. See the retained-GS guide.

`salmonpy.preparation.preflight_inputs(directory, manifest_sha256=None, registry=None)`
returns a JSON-ready input-input checking report without starting a solver. Explicit
trusted resolvers are optional; without one the scope is file integrity only.
See PREPARATION.en.md for batch-job exit codes and per-rank startup reports.

## Declared solver execution contracts

```python
from salmonpy.api import native_contract
contract = native_contract('fast_hse', route='shared')
assert contract['model_settings_owner'] == 'explicit model'
assert contract['runtime_verified'] is False
```

`native_contract` describes built-in `tddft`/`fast_hse` routes `single`, `shared`,
`registry`, and `fdtd` route `single`. It loads no library, executes no calculation
and enforces no new restriction. The result distinguishes single-state
ownership from shared/registered model ownership, caller-owned MPI, sequential main-thread
operations, entry-point collective ordering, current sign and failure scope.
Single-state standard-TDDFT C current is matter J=-j; its Python adapter returns charge
j. Fast-HSE and standard handle C currents are charge j. Follow exact headers for
query locality; a local diagnostic is not automatically collective.

This is declared documentation, not a live capability check or a stable-C-ABI
promise: `runtime_verified` is always false. Query actual runtime/checkpoint
capabilities separately. Unknown/custom component or route names return
`known=False`, `status='unverified'` and do not block a user-written calculation code; that
calculation implementation owns its contract. Each result is independent and may be annotated by
the calling program without changing later reports. Existing calls, numerical state,
current execution modes and accepted calculations are unchanged. Stable C/Fortran calling program
ABI conformance remains a separate milestone.

### Explicit live protocol probe

```python
from salmonpy.doctor import probe_native_protocols
report = probe_native_protocols('audited-build.json', component='tddft', route='shared')
# Inspect errors, observed protocols and protocol_match; runtime_ready stays unverified.
```

CLI: `python -m salmonpy.doctor --native-build audited-build.json --component tddft --route shared --json`.
Unlike `native_contract` and the separate static checker, this explicitly loads
solver code in a child process after checking the recorded library SHA256 (and
rechecking in the child). A hash proves identity, not trust or build-audit
compliance. Use only trusted, audited libraries: loading can run constructors.
The probe calls only known zero-argument version functions; it creates no solver
state and does not request MPI initialization/attachment, compilation, installation
or calculation. Unknown/custom routes do not even open the record or load code.

`protocol_match` describes the selected version queries only. Missing/newer
versions are reported, never converted to automatic calculation implementation rejection or feature
fallback. `runtime_ready` remains `unverified`, and operation signatures, physical
configurations and checkpoint/topology compatibility are not certified. A child failure
is reported without taking down the calling program; `library_loaded=None` means the failed
child could not report whether loading completed. Keep this opt-in diagnostic
separate from `check_job`, which never loads the supplied code.

Query-only C and ISO_C_BINDING Fortran hosts are in
`tests/native/test_native_protocol.c` and the `SALMONPY_ABI_PROTOCOL_PROBE` mode of
`tests/native/test_io_ownership.f90`. Development runner:
`tools/test_io_ownership.py BUILD_RECORD NEW_OUTPUT --abi-contract` (add
`--fast-hse` for a fast-HSE library). It checks runtime SHA, audits every linked
application object, compares returned versions and verifies MPI stays uninitialized.
This is version-query interoperability, not full model/state/step ABI conformance;
the C/Fortran calling interface remains experimental.

### Explicit state-lifecycle calling program test

`tools/test_io_ownership.py BUILD_RECORD NEW_BUILD_DIR --abi-state --fixture EXISTING_PZ_FIXTURE --state-output NEW_RESULTS_DIR --ranks 2` runs C and ISO_C_BINDING Fortran hosts against an audited standard-TDDFT library. Add `--fast-hse` with a HSE build record and an existing full-k HSE fixture to exercise the same checks. The HSE fixture explicitly sets k-only MPI decomposition for each calling program run; it does not modify the supplied GS. Add `--state-symmetry` for the reduced-k fixture: it references the supplied `sym.dat`, enables the symmetry-aware HSE initializer and drives z, the direction preserved by this fixture's retained operations. The runner checks that representative k points are actually reported. This is a fixture-specific direction, not a restriction on full-k x/y/z driving. The caller initializes/finalizes MPI. The hosts compare two live states with independent runs, reject a mismatched previous A without changing accepted state, reject stale handles after slot reuse, and confirm a fresh state starts independently. Existing GS files are read in place and checked unchanged; this test neither copies GS data nor writes electronic checkpoints. Select results storage explicitly.

The common directory-scope declarations are in `native_context.h`, included by both handle headers. On macOS, use the existing `tools/build_embedded_tddft.py::mpi_link_command` policy when linking MPI hosts too: private executable copies of OpenMPI Fortran COMMON sentinels can corrupt MPI_IN_PLACE operations across the shared-library boundary. The runner applies this policy and records link commands and executable hashes. Do not infer correctness from a query-only probe.

Verified scope: standard PZ and HSE full/reduced-k fixtures, C/Fortran hosts, MPI 1/2 ranks, paired versus independent currents. HSE uses its own documented signatures and rejection status (1 for mismatched previous A, versus 10 in standard TDDFT). The adapters live only in the common calling program test; they do not change either public ABI or numerical solver. The reduced fixture runs symmetry-aware PT-CN with learned exchange, 8 full k points represented by4 points. Other XC configurations, geometries and platforms require their own tests; the public C/Fortran calling interface remains experimental. These small fixtures check interoperability and state ownership, not converged material properties.

### Prepare the PZ calling program fixture

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


Additive API: salmonpy.preparation.prepare_fdtd(calculation, directory, mpi=None) writes a new FDTD called directly from Python input folder and returns its Path only. FDTDMPI.rank_arguments and command(executable) configure/preview argv, without execution. See PREPARATION.en.md.

Preparation-only API: prepare_dft(model, directory, scf=None, save_restart=True) and prepare_tddft(input_file, directory, files=None) return a Path without SCF, TDDFT workers, MPI initialization or submission. See PREPARATION.en.md.


## Review pause for unconverged SCF (2026-11-21)

`DFT.run()` and `NativeGSSession.run()` warn and raise `SCFReviewRequired`
when SCF exhausts its iteration cap without convergence. Automatic downstream
computation stops for review. Completed data are retained: `exception.result`
is the saved batch result or still-live calculation object; `exception.report_directory`
points to residual CSV, optional PNG, and report.json. Missing Matplotlib does
not prevent CSV output. A report-writing failure still pauses continuation and
warns to inspect the original solver output. Normal GS saving settings remain;
`save_restart=False` does not collect a new GS restart dataset.

```python
from salmonpy.errors import SCFReviewRequired
try:
    gs = model.run(engine, 'new-gs-output')
except SCFReviewRequired as review:
    gs = review.result
    print('Review SCF residuals:', review.report_directory)
    # End here and inspect residual trends, threshold and input conditions.
    # Do not automatically start TDDFT from this handler.
```

Inspect decreasing, stalled or oscillating residuals before adjusting iteration
or mixing settings. Only explicit user approval after review warrants
`allow_unconverged=True`; an AI must not automatically retry with this option.
`GroundStateResult.open()`, reading/saving diagnostics, and warm starts remain
available. `band_initial_state()` and `tddft_initial_state()` also pause for a
known unconverged batch GS unless the reviewed override is supplied.

Intentional solver step budgets (`run(..., steps=K)` and
`continue_scf(steps=M)`) returning `paused=True` are normal partial SCF states.
Only terminal nonconvergence pauses for review. The calculation object remains live
for diagnosis or saving. All solver ranks must agree on `allow_unconverged`.
Numerical update order and SCF thresholds are unchanged. This replaces the old
policy of automatically proceeding after an unconverged-SCF warning.


## FDTD/multiscale energy budget

`analysis.absorbed_energy_density(time, E, j, units='atomic_gaussian',
current_convention='charge', field_kind='total')` integrates collocated total
E dot physical charge current j per macro point in Fortran (trapezoidal rule).
E/j have shape (time, point, 3); output is Hartree/bohr³. Explicitly convert
SALMON matter current J to j=-J first. This is net work including reversible
material storage, not heat alone. Do not add LD pole storage or electronic
energy change again when counting this work.

`analysis.reflection_transmission` takes separated incident/reflected/transmitted
E/H arrays, each (time,3), and returns R(omega), T(omega) and a validity mask.
The normal points from incident to transmitted side. Ratios use complex Poynting
flux, retaining polarization and impedance, not an E-only amplitude ratio.
Without area_weights, scope is uniform equal-area plane ports. Separate wave fields with a reference run
or a directional decomposition first; a total field is not a reflected field.
Collocate Yee E/H in both space and time. Weak-drive bins remain invalid/NaN.

```python
from salmonpy.analysis import absorbed_energy_density, reflection_transmission
from salmonpy.fdtd import energy_balance

work = absorbed_energy_density(time, E_points, charge_j_points,
    units='atomic_gaussian', current_convention='charge', field_kind='total')
ports = reflection_transmission(time, E_inc, H_inc, E_ref, H_ref, E_tr, H_tr,
    normal=(0, 0, 1), units='atomic_gaussian', window='none')
check = energy_balance(work, point_volumes, ports,
    material_ids=material_ids, transverse_area=area,
    field_energy_change=delta_field_energy, other_outgoing_energy=other_outflow,
    units='atomic_gaussian', rtol=.01)
print(check.residual, check.relative_residual, check.within_tolerance)
# ports.omega/R/T/valid; work retains the macro-point order.
```

Residual = incident-reflected-transmitted-sum(work*volume)-field_change-other_outflow.
Volumes are bohr³, area bohr², total energies Hartree. Preserve material IDs;
repeated IDs still refer to independent macro points. Required field_change
includes instantaneous epsilon_infinity field storage, excluding material work
already counted. Include unmeasured lateral flux/PML absorption separately,
with one consistent region boundary; never count the same outgoing energy twice.

Spectra use the selected window; total-energy fluences always integrate raw,
unwindowed E/H. Do not average R/T or infer all absorption from drive-band
1-R-T alone. Time-domain flux also includes frequency conversion. Truncated
traces, residual fields and delayed emission can affect the budget. The default
1% is advisory, not accuracy certification. Excess warns without stopping or
changing a solver. Zero incident energy and nonfinite diagnostic inputs are
rejected. After restart, join traces without duplicate boundary samples and
compare the same initial/final states. No added MPI communication or implicit
solver execution. Nonuniform periodic planes use the area_weights entry below; arbitrary closed surfaces are not automatically collected.

Stored-array example: `python -m salmonpy.examples.fdtd_energy INPUT.npz NEW_OUTPUT`.
See its module docstring for NPZ keys/units. It never launches a solver.


### Real coupled-run checks (2026-11-21)

On Apple Silicon with an audited SALMON executable, absolute budget residuals
relative to incident energy were 0.0108% for solver LD+E/H FDTD, 0.00270% for
an already-running half-spacing/half-dt LD check, and 0.000711% for solver
TDDFT+Weyl FDTD. These verify coupling/observation accounting, not converged
real-material optical properties. No further precision search was performed
after the user's acceptance of this accuracy.

The LD fixture is a transversely periodic uniform slab. Current was reconstructed
from measured solver E with the identical solver LD calculation routine; final current matches
solver LD checkpoint data at rtol=1e-10, atol=1e-18. Raw Yee H was centered in
space; the saved time-centered H was not averaged again in time. TDDFT used
an existing GS in place, 2 macro points, 8 MPI ranks and dt=.003 au. E/H/j and
field energy were read directly from solver output files; no new GS/wavefunction data
were saved.

The primary check delays incidence so the initial field is empty. A separate
prefilled-wave check includes the initial field energy in the budget. Finite-time
R/T can appear above one when incident traces are truncated; inspect initial
fields and observation coverage rather than clipping the spectrum.

Source-checkout reanalysis tool: `tools/analyze_fdtd_energy.py`. It only reads
completed validation fixtures and explicitly builds/audits a small C reference;
it does not launch SALMON or MPI. Wheel users can run `examples.fdtd_energy`
on the exported NPZ observations. Real-material runs are not added to required CI.


### Complete transverse-periodic planes, including diffraction

Pass `area_weights=dA` to `reflection_transmission` for nonuniform planes.
E/H have shape (time, Ny, Nz, 3); dA has shape (Ny,Nz) in bohr². Flattened
(time, patches,3) arrays also work. Cover one complete transverse periodic cell
and avoid double-counting a repeated endpoint. Compute local E cross H BEFORE
summing areas: crossing mean E and mean H loses spatially phased diffraction.

```python
from salmonpy.analysis import surface_flux, reflection_transmission
from salmonpy.fdtd import energy_balance

# Each field: (N, Ny, Nz, 3), already collocated in space and time.
# dA: (Ny,Nz); each entry dy*dz for a uniform grid.
ports = reflection_transmission(time, E_inc, H_inc, E_ref, H_ref, E_tr, H_tr,
    normal=(1, 0, 0), area_weights=dA, units='atomic_gaussian')
area = ports.config['surface_area']
check = energy_balance(point_work, point_volumes, ports,
    material_ids=material_ids, transverse_area=area,
    field_energy_change=delta_field_energy, units='atomic_gaussian')
net = surface_flux(time, E_total, H_total,
    normal=(1, 0, 0), area_weights=dA, units='atomic_gaussian')
# net.power: signed total Hartree/atomic_time; net.energy: cumulative Hartree.
# net.fluence: area-averaged Hartree/bohr².
```

R/T likewise sum area-weighted complex Poynting flux before taking ratios.
All spatially resolved orders contribute, including oblique waves and modes
whose mean field vanishes. Subtract reference E AND H to form a reflected field;
never substitute an intensity difference. For oblique diffraction, use measured H,
not normal cross E. Put planes clear of material/PML, with matching cell areas.
Finite-time residual/near fields and weak drive can give unstable signed ratios:
inspect masks, warnings, duration and plane positions; do not clip negative bins.

With matched transverse periodicity, opposite lateral fluxes cancel over one
cell. Front/back integrated flux, material work and field-storage change close
the budget. Extra lateral sources/nonperiodic channels must be counted separately.
No added MPI calls: disjoint patch power/energy add on the owning group. Do not
average rank-local R/T; evaluate ratios over the whole cross section. Automatic
closed-surface placement, order-resolved decomposition and spatial interpolation
are not provided. Batch FDTD can now collect planes through PlaneMonitor below. External collocated
arrays remain supported; the old point-observation entry retains its200-point limit.

Fortran performs area/time integration. The legacy uniform-port path is unchanged.
The budget checks transverse_area against the recorded plane area to prevent
area double-counting. The stored-NPZ fdtd_energy example accepts area_weights.

A real solver transversely periodic lossless dielectric grating (z invariant)
closed at about0.0107% of incident energy with32 transverse samples. Averaging
fields first lost transmitted energy equivalent to about10.1% of incidence.
Final stored field energy was included from solver field checkpoint data. This
is local observation/accounting verification, not general3D material accuracy
or multi-node scaling certification.


### solver acquisition and automatic E/H colocation

Specify `planes=(PlaneMonitor(axis, position), ...)` in batch FDTD. Axes are0/1/2,
positions bohr. The nearest grid plane is selected and requested/actual positions
are recorded. Transverse faces must be periodic; keep the plane and preceding
neighbor clear of normal boundaries/PML. This stream is separate from the old
200-point observation entry and does not expand a large plane into point probes.

```python
from dataclasses import replace
from salmonpy.fdtd import PlaneMonitor
from salmonpy.geometry import Geometry

# calculation already has transverse periodic boundaries. Construct/prepare only
# on a management node; run explicitly in the permitted compute environment.
monitored = replace(calculation,
    planes=(PlaneMonitor(0, x_front), PlaneMonitor(0, x_back)))
reference = engine.run(replace(monitored,
    geometry=Geometry(monitored.geometry.grid), materials=(), smoothing=None),
    'new-vacuum-reference')
result = engine.run(monitored, 'new-material-run')
plane = result.planes[0]
print(plane.time, plane.E.shape, plane.position)
net = plane.flux()
ports = result.reflection_transmission(reference)
check = result.energy_balance(reference, point_work, point_volumes,
    material_ids=material_ids, field_energy_change=delta_field_energy)
```

PlaneData E/H have shape(time,N1,N2,3), area_weights(N1,N2). Fortran averages
2 neighboring E edges and4 H faces at a common position with periodic wrapping.
solver saved H is already centered at E time: no second temporal average.
Only ranks owning the plane/preceding plane append per-plane streams. Offline
tile assembly and colocation add no MPI collectives. This first entry is fresh
batch FDTD, with explicit restart-segment reading/joining as described below;
in-process NativeEHField uses explicit sample_planes snapshots, described below;
existing numerical checkpoint/restart contracts are unchanged.

Times include0 and each observation stride. If steps is not divisible by stride,
the last observation is the last stride boundary. Match material work and field
storage to this interval. They cannot be inferred as zero from plane data and
remain explicit inputs. Use matched reference grid/source/clocks/plane placement.
Returned fields are immutable. `read_native_planes(directory, calculation,
ranks=...)` rereads completed streams without a solver or MPI.

A matching runtime contains SALMONPY_NATIVE_PLANES_V1. Old runtimes or failed/
incomplete plane output warn and mark diagnostics unavailable while preserving
a completed FDTD calculation and its point outputs. The budget never fabricates
missing planes. Imports/preparation do not build; explicitly obtain/build a
matching runtime if needed. Raw field output is approximately
8*6*2*N1*N2*(1+steps//stride) bytes per plane, plus tags/tile headers.

Actual grating tests preserved all128 prior point-output files and full EM
checkpoint bytes with monitoring off/on. Two transverse MPI ranks gave matching
automatically acquired planes/budgets (~0.0107% residual). A separate normal split
placed the plane/preceding plane on different ranks and retained bit-identical
nonzero colocated E/H. These are local implementation checks, not multi-node
scaling or converged real-material optical accuracy.

### solver batch restart observation segments

The existing solver restart loader/solver is unchanged. A restarted stream records
only samples strictly after the checkpoint, at the original global stride; do not
insert a zero sample or average H again. solver batch restart keeps the original
`nt_em` and observation stride. Launch the restart explicitly on the compute node
using existing SALMON input files; these reading functions never launch a solver.

```python
from salmonpy.fdtd import read_native_planes, join_plane_segments, FDTDResult
from pathlib import Path
# calculation retains the original total steps and geometry.
# early_dir contains a complete stream ending at the checkpoint (not a full run).
early = read_native_planes(early_dir, calculation, ranks=2, stop_step=checkpoint_step)
late = read_native_planes(restart_dir, calculation, ranks=2, start_step=checkpoint_step)
planes = join_plane_segments(early, late)
result = FDTDResult(Path(restart_dir), (), planes)
budget = result.energy_balance(reference, work, volumes,
    material_ids=ids, field_energy_change=stored_delta)
```

`stop_step` defaults to the original final step. All expected samples through
that step must exist and trailing data are rejected: archive a complete segment,
never pass a still-growing file. Joining verifies plane geometry, global clock,
stride, total planned steps and exact checkpoint continuity; gaps, overlaps and
configuration differences are errors in the diagnostic only. No replay, sample
invention or implicit removal of duplicates occurs. Source arrays stay unchanged
and joined arrays are read-only. Reference/work/storage must cover the joined
interval. in-process sampling/publication uses the explicit committed-only entry below.

### in-process FDTD: committed-only plane snapshots

Set the same `FDTD.planes`, then all field-group ranks explicitly call
`field.sample_planes()` at step0 or an `observation_stride` boundary. The solver
`sp_eh_plane_read` ABI reads only owned saved E/H tiles, without MPI, file writes
or mutation. Python performs one field-group gather containing tiles and errors;
only field rank0 receives immutable PlaneData. Other ranks return None. No world
collective, per-step automatic communication or replay is added. H is already
centered at E time. An older library or unavailable diagnostic warns on field
rank0 and returns None, leaving the calculation usable.

Call after `step_protected`/protected restore has committed and released its
candidate. During `_step_guard`/restore-candidate ownership, sampling/publication
warns and returns None, without exposing a speculative frame or creating files.
`output_directory=NEW_DIR` optionally saves that single snapshot via the existing
no-overwrite output archive. Archive failure warns and still returns samples;
existing archives and field state are preserved. This is a separate immutable
observation, not a live stream automatically rewound by future restores. For a
new branch use a new destination. Crash/OOM/MPI failure atomicity is not claimed.

```python
from salmonpy.fdtd import combine_plane_samples, FDTDResult
# Every field rank calls, only at committed steps and solver observation times.
snapshots = []
for step in range(steps):
    coupled.step_protected(targets)  # called by all participating groups
    if field is not None and (step + 1) % snapshot_interval == 0:
        sample = field.sample_planes()
        if sample is not None:
            snapshots.append(sample)  # field rank0 only
if snapshots:
    planes = combine_plane_samples(*snapshots)
    result = FDTDResult(output_directory, (), planes)
    budget = result.energy_balance(reference, work, volumes,
        material_ids=ids, field_energy_change=stored_delta)
```

Snapshots must form a matching, strictly increasing uniform time series; no gap
filling/duplicate dropping occurs. Use at least two samples for flux/spectra, a
sample interval resolving the fields and a reference on identical times. Align
explicit material work/storage to this observation interval. Storing snapshots
costs memory proportional to plane area × number of samples; select a stride.
`mpi_mixed_native --planes-x FRONT BACK --plane-stride 16` explicitly selects
transverse-periodic boundaries, captures only after completed steps, and writes
one `plane-series` archive at the end. Existing default sample behavior is unchanged.

Local LD fixture: direct/protected restore/retry and serial/MPI2 planes preserve
traces and point-output bytes. Across the prior batch build E/H agree within
1e-12 relative/1e-14 absolute (maximum observed difference below8e-19). Numerical
checkpoint payloads match; legacy unused coordinate padding may differ across
independent processes, so full-file equality is not claimed for every MPI rank.

Coordinate-persistence update (2026-11-21): new runtimes zero unused coordinate
padding before unchanged coordinate/halo arithmetic. Serial/MPI2 currents, fields
and point outputs match the prior build; complete checkpoint files now match on
every rank after protected/direct restore. The padding difference above records
the pre-fix validation, not an unresolved issue in the new runtime.


## Common operations: arguments, lifetime and failure

`native_contract(component, route=...)['host_operations']` describes the six
common operations below. The same declaration sets Python's C argument/return
types. All required symbols are resolved before any operation executes; a missing
symbol on one solver rank rejects the group before attachment. This uses the
existing admission agreement and adds no timestep communication. Custom backends
keep their own contracts. Describing a signature does not inspect machine code
or certify every operation of a solver. Use authenticated builds and host tests.

| C entry | Arguments → result | Communication / precondition | Effect and failure |
|---|---|---|---|
| `sp_attach` | `int64_t MPI_Fint` → `int` | collective on supplied live intracommunicator; no attachment | duplicates communicator; 0 success, −1 handle width/range, otherwise MPI status |
| `sp_detach` | none → `int` | solver-group collective; close all owned states first | frees duplicates; never finalizes host MPI; 0 or MPI status |
| `sp_info` | `int* rank, int* size` → void | local; attached | writes caller buffers; no returned error status |
| `sp_barrier` | none → void | solver-group collective; attached | same call order on every rank; no returned error status |
| `sp_sum` | `double` → `double` | solver-group collective; attached | sum returned to every group member; no returned error status |
| `sp_library` | `char* buffer, int capacity, int* active` → void | local; allowed before attachment | capacity ≥1; truncated NUL-terminated MPI text; active=1 only while MPI is initialized and not finalized |

Every call is sequential on the main thread. All pointers must be valid and
owned by the caller; buffers are not retained after return. Keep the supplied
communicator live through attachment. Info/sum/barrier do not change field or
electronic state. Arbitrary MPI/process failures and invalid pointers are outside
the rollback guarantee; void/value-returning entries cannot report recoverable
errors. The Python context rejects operations after close and repeated close is
a Python no-op; raw C detach still requires a live attachment.

The current solver-specific operation declarations remain in
`standard_tddft_handles.h`, `standard_tddft_models.h` and `fast_hse_handles.h`:
creation/register/spawn, step/info, checkpoint/transaction and close have their
own collective order, token, current-sign and profile requirements. The common
operation table does not replace or certify those contracts. Full operation and
platform coverage remains on the task list.


## TDDFT and fast-HSE state operations

For the built-in current interfaces,
`native_contract(component, route=...)['state_operations']` lists required C
symbols, argument order, result, collective scope, ownership and mode meanings.
The same declaration supplies Python bindings for single/shared/registry routes.
Resolve every required symbol before binding or state execution; retain separate
protocol/profile admission. Optional timers, electron-number diagnostics and
retained-GS imports keep their own capability queries. They are not implied by
this required-operation list. Custom backends remain independent.

Arguments are C values unless marked pointer. `int64` is an opaque signed
64-bit token, not an address; `string` is a NUL-terminated filesystem byte string.
`vector3` means contiguous three-element float64 Cartesian data,
`decomposition5` contiguous five-element int32, `history8` eight-element int32.
`int*`, `int64*`, `double*` are scalar output pointers; `buffer` is writable
caller storage. All buffers are caller-owned and not retained after return.
Protocol queries are local; state calls are main-thread, sequential collectives
on the owning solver group in the same order, with the model's I/O directory.
Only successful outputs may be consumed. Nonzero status does not universally
mean rollback succeeded; use the operation-specific recovery contract.

| Operation | Lifetime / behavior |
|---|---|
| handle_init / model_register | Return a positive owning generation-checked state token; registration also returns a negative model token. Keep inputs/GS/PP immutable. Close each state exactly once. |
| model_new / model_spawn | Start an independent state from original GS, never from a sibling excitation. TDDFT inherits dt/capacity/decomposition/A0; HSE receives dt/capacity per call. |
| handle_step | Old and new A/c endpoints in atomic units; old must exactly equal the accepted endpoint. Standard handle APIs and HSE return charge j; standard single-state C API returns matter J=−j. HSE also writes total energy in Hartree. |
| handle_info / handle_label | Query accepted step/A/j or a state-owned output prefix. Label capacity includes NUL; too small returns13. |
| handle_checkpoint | Modes0 save,1 fresh load,2 stage,3 commit,4 cancel. Paths may differ by rank. Caller authenticates runtime/input/topology compatibility. Raw load alone is not whole-session atomic restore. |
| handle_transaction | Standard:0 capture,1 rollback/close streams,2 release,3 reopen. HSE:0 capture,1 rollback/close streams,2 reopen/flush/release (retain on failure),3 reopen. Caller restores only this owner's stream bytes between rollback and reopen. |
| handle_close / model_close | Consume state token / retire model allocator. Existing children retain model resources until last close. A stale token stays invalid after slot reuse. Close owned states before detach. |

Shared/registry rejection codes include11 mixed rank arguments,12 stale token
or attachment/lifetime mismatch,13 buffer/capacity limits; other codes depend on
operation/profile. See the shipped C headers for exact prototypes. The required
single route list covers the existing initialization/update/shutdown entry
points; it is not a declaration that all single-state optional operations are
covered. No universal recovery from MPI failure, STOP or invalid pointers is
promised. This work standardizes current signatures; broader profiles/platforms
still require separate evidence.

Each entry also records ordered `parameters` with name/type/in-or-out.
`decomposition_order` is k/orbital/grid-x/grid-y/grid-z. Standard single-state
C output is explicitly named `J_matter`, distinguishing it from charge j.

## FDTD field operations

`native_contract('fdtd', route='single')['field_operations']` describes current
FDTD arguments/directions, layout, communication scope and restore/transaction
modes. Python uses the same declarations; C prototypes ship in `native_fdtd.h`.
Only init/step/close are required at field creation. Ownership routing, H and plane
sampling, checkpoint and protected-step entries are checked when used; they do
not restrict custom backends or unused optional features.

Field arrays are C-contiguous float64 `(points,3)` (Fortran `(3,points)`). Init
indices are int32 one-based global triples; Python converts zero-based indices.
The C argument named J contains physical charge-current density j at the previous
half-step, in atomic units; no sign reversal or cell-volume multiplier. Checked
steps also receive accepted A/c. Owned steps pack J in ascending current-mask
row order; A/E still contain every mapped point. Fields, LD and update order are
unchanged. Use COUPLED's E/H/A time-placement contract.

All state calls are sequential on the main thread. Init/step/close, H sampling,
checkpoint and transactions are field-group collectives. `current_mask` is local,
and `plane_read` is a local read of committed sampled data: axis1..3, local-grid
offset, int32 info[21], first call with count0 for metadata, then float64 values.
Protocol queries are local and do not initialize MPI or the field.

Checkpoint modes0/1/2/3/4 mean save/load/stage/commit/cancel, with distinct paths
per rank. The caller validates build/input/topology and coordinates publication;
raw save is not a multi-file atomic publish. Field transaction modes0/1/2/3 mean
capture/rollback/release/no-op, unlike electronic stream-reopening transactions.
Unsupported observation/output profiles return15 for protected-step capture.
Buffers remain caller-owned; consume outputs only on success. No guarantee covers
arbitrary process/MPI/OOM failure.

## GS and SCF operation contracts

`native_contract('gs',route='single')['gs_operations']` describes required and
optional GS entries with ordered names/types/directions, ownership and modes.
Python bindings use the same declaration. `single` selects the GS interface,
not a one-state limit: every live GS/snapshot has its own owning token. Close
snapshots/owners before detaching. Pure protocol/format and resource-count queries
are local; state info and array reads are collective because tokens are agreed.

Create runs SCF; start/continue preserve initialization and true mixing history.
An unconverged final iterate is retained, but Python pauses automatic continuation
for residual review by default. `allow_unconverged=True` is an explicit reviewed
choice; the finite 8-iteration developer fixture has its own explicit CLI flag.
It does not demonstrate physical convergence. Import reads existing GS without
SCF; unavailable convergence/forces must not be invented. Bootstrap for direct
SCF restart initializes static resources while skipping numerical initialization
and SCF iterations.

array_info returns ndim and int64 shape/lower[7]; read buffers are float64 in
Fortran order, can include local halos, and must match the queried extent. info
uses int32 stats[4], float64 values[3] (Hartree energy/residual/tolerance);
scf_info uses int32[6] including lifetime initialization/execution counters.
In-memory snapshots belong to the same live owner: restore copies numerical/XC/
mixing/HSE history without replay, but does not rewind text logs or lifetime
counters and is not a portable disk archive.

Candidate modes2/3/4 stage/publish/discard; publish/discard use snapshot0. The
host must agree before publication. Disk mode0 writes a new rank-local file,
mode2 stages a read; candidate(state,0,3) publishes only after authentication.
Repartition entries require their own layout-specific capability versions and
ordered old-rank bounds. No unused feature is required and no new hard stop is
added for valid/custom calculations. Optional timers/GS-to-RT bridges still use
separate queries. `doctor --component gs --route single` checks only the common
host descriptor and GS protocol, without starting SCF or certifying all operations.

## Unit conventions and NLCC validation

`unit_contract()` declares the current conventions without loading a solver:
boundary quantities are atomic units; SALMON's retained eV input constant is
27.21138505 eV/Hartree, geometry uses0.52917721067 angstrom/bohr, while Python
energy presentation uses27.211386245988 eV/Hartree. The energy constants differ
by4.3951750e-8 relatively (1.195988e-6 eV per Hartree). ASE forces retain the
existing presentation-energy/geometry-length ratio. Do not silently replace
input constants with display constants: that would alter interpreted inputs.
The separate legacy +U constant27.2116 remains unchanged and is outside the
current declared GS/TDDFT profile. This contract documents current behavior; it
does not inspect or certify a custom solver's unit conventions.

NLCC forces require binary-bound evidence from the maintained source tree. The
current tree hash must match both runtime and receipt/native evidence; missing
source identity or old receipts reject without fallback. Current prescribed
tests cover6 core-force geometries and both MPI core-table formats, including
2/4-rank ×1/2-thread replication and collective nonfinite rejection. An MPI-linked
validated executable additionally requires matching MPI evidence. This is
implementation validation, not all-material or GPU physical certification.

The maintainer source-checkout tool `tools/validate_nlcc_runtime.py` accepts `--mpi-record` and an explicit completed
`--application-record`, avoiding repeat DFT when the whole-application scalar
comparison already exists for the exact binary/reference. Evidence alone is
not publisher authentication: use trusted build/validation records. Uncertified
executables do not gain NLCC forces. Current validation used the fixed Si PBE
reference: energy difference0, maximum force difference9.773e-16 Hartree/bohr.
No numerical formula/constants or reference values were changed.

## AI task discovery

`salmonpy.manual.search(query, lang="en", limit=5)` returns ranked summaries without guide bodies. `task_details(task_id, lang="en")` returns one exact task with its requested-language guide. If unavailable, it returns the source guide with `translation_available=False`. `load()` reads the entire JSON; prefer search for routine guidance.
