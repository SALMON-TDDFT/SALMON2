# Full solver runtime

LD/SBE wheels already include their compiled library. DFT, TDDFT, FDTD and band
calculations use a **linked executable plus its runtime libraries**. An object
file alone would still require a linker and dependencies on the user's machine.
Importing SALMONpy never downloads or builds anything.

## Compute nodes without internet access

Do not assume compute nodes can reach external services. Separate acquisition
from execution: prepare files on an internet-connected workstation or a preparation
host permitted by the facility, then transfer them. Login nodes may also be offline;
follow the facility's access policy.

1. **Identify the target first:** record OS, CPU, Python version, required features,
   MPI implementation and site modules. A wheel/runtime for a local Mac will not
   run on a Linux compute node.
2. **Stage all dependencies:** obtain compatible SALMONpy/NumPy wheels and optional
   dependencies such as ASE, required runtime ZIPs and trusted catalogs, structure/
   supercell files, pseudopotentials, inputs and initial states. Include referenced
   files; record versions, sources, licenses and SHA256 hashes. Check which MPI,
   BLAS, Libxc or other libraries are bundled and which the site must supply.
   Do not embed access credentials in calculation inputs.
3. **Verify after transfer:** check hashes and preserve relative references. Install
   Python packages from the local wheel collection with package-index access disabled.
   Install runtimes with both a local `--catalog` and `--archive`, then run `verify`
   and a small example on the actual compute node. Success on a login node alone
   does not establish compute-node compatibility.
4. **Resolve missing files before submission:** use downloaded structure files, not
   live database URLs during execution. Report missing dependencies explicitly;
   do not silently download or switch to source builds. If no compatible runtime
   exists, explicitly build in a site-compatible environment with dependencies
   already staged.

This guidance is not an automatic offline bundle/dependency resolver. Custom
backends must also be prepared to execute without network acquisition. SALMONpy
import is offline; it does not sandbox network calls made by user callbacks.

## Prefer an existing or prebuilt runtime

First establish the user's memory limit, CPU budget and desired workflow. Check
OS/architecture and required features; LD wheel compatibility does not establish
full-runtime compatibility. Reuse a previously audited executable via
`SalmonExecutable(path, recorded_sha256)`, or verify a managed runtime:

```python
from salmonpy.runtime import load_runtime
engine = load_runtime('my-runtime', required=['tddft', 'libxc'])
# Pass executable=engine to the existing DFT/TDDFT/band loaders.
# For FDTD: engine.run(calculation, 'new-output', threads=2).
```

For public releases, use the trusted publisher's catalog URL:

```sh
python -m salmonpy.runtime status --catalog CATALOG_URL --require tddft libxc
python -m salmonpy.runtime install my-runtime --catalog CATALOG_URL --require tddft libxc
python -m salmonpy.runtime verify my-runtime --require tddft libxc
```

The [independent SALMONpy preview Release](https://github.com/otobetomohito/SALMONpy/releases/tag/runtime-2026.11.14) is available and the bundled catalog registers three platforms. This repository is private: download `catalog.json` and the matching OS/CPU runtime ZIP through an authenticated browser, then use `--archive` below. Anonymous URL downloads cannot access private assets. This Release does not include Libxc. `CATALOG_URL` above illustrates a public publisher only.
Multiple matches require `--artifact ID`. Destinations must be new. A missing
match raises `RuntimeUnavailable` and prints the explicit build entry point;
it never silently drops Libxc, changes architecture or starts compilation.
Downloads use HTTPS and SHA256. Archive paths and all installed files are checked.
Catalogs must come from a trusted publisher: their hashes are not signatures.

## Current branch artifacts

The **Full runtime bundles** Actions workflow builds macOS arm64, macOS x86_64
and Linux x86_64 profiles and uploads a bundle only after relocation and native
TDDFT tests succeed. A configured runner is not evidence of a passing run. These
profiles use a single MPI rank and PZ/TB-mBJ; they do not include Libxc. Windows
full-runtime bundles are not yet provided (Windows LD/SBE wheels are separate).
The locally tested macOS arm64 Libxc bundle is a separate artifact/profile.

Download the matching successful run's `salmonpy-runtime-...` artifact from the
current branch's Actions page. Extract the **outer Actions ZIP**, retaining the
inner runtime ZIP and `catalog.json`. Then:

```sh
python -m salmonpy.runtime status --catalog downloaded/catalog.json --require tddft
python -m salmonpy.runtime install my-runtime --catalog downloaded/catalog.json --archive downloaded/RUNTIME_ZIP --require tddft
```

`RUNTIME_ZIP` means the actual ZIP filename inside that artifact. CI catalogs
have no public asset URL; installation requires this local archive. No compiler
or linker is needed. OS versions older than the recorded build baseline are
rejected conservatively. Native calculations still need compatible inputs and
pseudopotentials, which are not part of the executable bundle.

## Explicit fallback build

Only offer a build when no appropriate prebuilt/existing runtime is available.
Explain required tools and confirm the user's CPU/memory budget. Use the bundled sources and explicitly build in a new directory:

```sh
python -m salmonpy.runtime build NEW_BUILD --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=gfortran \
  --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DBLA_VENDOR=OpenBLAS \
  --cmake-arg=-DOPENMP_FLAGS=-fopenmp \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --cmake-arg='-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -ffree-line-length-none -fallow-argument-mismatch'
```

These are GNU compiler options, not universal platform flags. Supply the actual
GNU C compiler path on macOS; `/usr/bin/gcc` is normally Apple Clang. Install
CMake >=3.22, GNU Fortran/C and LP64 OpenBLAS beforehand. The command never installs
system dependencies. Windows native full-solver support is unvalidated; use a
supported Linux environment/WSL or a separately validated site build.
`--jobs` controls build concurrency, not simulation threads or a memory limit.
For Libxc supply `--cmake-arg=-DUSE_LIBXC=ON` and
`--cmake-arg=-DLIBXC_INSTALLDIR=LIBXC_PREFIX` for a compatible Fortran-enabled Libxc.

The wheel includes the complete pinned solver sources, verified file inventories and audit scripts,
so no SALMONpy checkout is needed. Source preparation, compilation and object
checks are logged in `NEW_BUILD/build.log`; `build.json` records status, exact
commands and the resulting executable SHA256. Failures preserve the log and
source; retry in a new directory. Success means **built and object-audited**,
not independently validated scientific accuracy. Run the relevant native
regression before production and use the recorded hash:

```python
import json
from pathlib import Path
from salmonpy.fdtd import SalmonExecutable
record = json.loads(Path('NEW_BUILD/build.json').read_text())
if record['status'] != 'built-and-audited':
    raise RuntimeError('runtime build did not complete')
engine = SalmonExecutable(record['executable'], record['sha256'])
```

MPI runtimes must match the cluster's MPI stack and remain explicit site builds.
DFT/bands and process-isolated TDDFT retain their single-rank scope. FDTD can use its existing spatial MPI decomposition through FDTDMPI; see FDTD_MATERIALS.en.md. This acquisition
layer does not alter the DFT/TDDFT update order, state isolation or A→J interface.

Experimental host-MPI TDDFT and native E/H coupling now have a separate entry: [COUPLED.en.md](COUPLED.en.md). One point per solver group and one native lifetime per process; do not infer unrestricted reentrancy or restart support.

## Validated NLCC force capability

`load_runtime(directory, required=['nlcc-forces'])` enables Calculator NLCC forces
only with checksum-verified evidence matching the binary, supported patch and
CPU/unpolarized/NC scope. A feature name or ordinary executable constructor alone
is insufficient. For a trusted local build, use
`load_validated_build(build_json, validation_json)`. See [Calculator](ASE_CALCULATOR.en.md).

Developers explicitly run XC/core tests with `tools/run_nlcc_native_tests.py`, then
`tools/validate_nlcc_runtime.py BUILD_JSON XC_TEST_JSON CORE_TEST_JSON PBE_UPF NEW_CHECKS_DIRECTORY`.
For a bundle, pass `--nlcc-validation CHECKS/validation.json --nlcc-potential PBE_UPF`
to `tools/package_salmon_runtime.py`. The relocated executable reruns the application
regression and receives its own binary binding. Import never builds or runs a solver.
This is implementation regression evidence; material grid/k-point/SCF/force
convergence remains the caller's responsibility.

## Experimental RT-owned step configuration

The default runtime is unchanged. Explicitly apply `_runtime_source/tools/prepare_rt_step_config.py BASE NEW_DIRECTORY` to a matching externally prepared runtime; in a development checkout use `python tools/prepare_rt_step_config.py BASE NEW_DIRECTORY`. Only a new external directory is created. Inputs/fixed executables remain unchanged; failed preparations are retained with status and are not approved for execution.

This first extraction gives each RT state its dt, Taylor order, propagation and field-integration settings. Validation uses synthetic PZ inputs; full native multi-material contexts are not implemented. Isolated calcVbox, ionic motion, single-scale Maxwell and other shared settings remain outside this milestone. FDTD coupling still requires a common dt.

The user's fast PR-intended implementation is the HSE baseline. This PZ regression is not an HSE speed/accuracy benchmark and the overlay is not automatically applicable to that PR branch. HSE ACE factors, prediction histories and mutable workspaces need separate ownership/validation. Explicit builds must resolve matching dependencies/compiler flags and retain source/all-application-object prohibition audits and CMake checks.

## Experimental fast-PR HSE RT state extraction

A separate HSE profile pins the user's PR-intended revision `c9e8ce431f43eefd4a88998134adb4796f8ec837`. Run `python tools/prepare_hse_rt_state.py NEW_DIRECTORY` to prepare a new external source copy, or add `--reference-only` for the unchanged numerical reference. Both retain mandatory source/object audits in CMake configure/incremental builds. The PR checkout and default SALMONpy runtime are unchanged.

This stage groups ACE factors, accepted prediction histories, PT-CN history/dt, Taylor trials and action caches under `s_exx_rt_state`. Lifecycle operations take explicit state arguments and reject invalid dt/stages before mutation. Optional owned-allocation byte accounting includes allocated history-array descriptors/fixed fields and nested buffers; it is not RSS.

The legacy HSE adapter still uses one default state. Geometry/kernel/Wannier resources, public legacy flags and input settings remain shared, so this does not provide full native multi-material contexts or concurrent calls. Fast exchange paths are retained; synthetic numerical regression does not establish real-material accuracy or a speed improvement. FDTD coupling still requires a common dt.

## Experimental explicit HSE RT routing

Add `--explicit-routing` to the preceding preparation command to copy the completed routing sources after verifying their input/output hashes. It cannot be combined with `--reference-only`. Default preparation and existing calls remain unchanged.

Twenty-four exchange/history procedures forward a caller-supplied RT state through their internal calls. Resolution uses local references, without switching a process-wide active state or copying all arrays. Public legacy flags/energy are state-owned; old names are aliases to the default owner only. Reassociating those aliases is unsupported.

Without a system binding, the state-only routing scope remains common exchange configuration, serial CPU scalar dense ACE action and lifecycle. Nondefault refresh, direct/full, packed/distributed and spin routes are rejected before shared resources are reached. With optional `rt_status`, `exx_refresh`/`exx_add_action` return nonzero for this capability mismatch without changing accepted state/output. Without it native errors may terminate the process. This is not general failure atomicity or material C ABI conformance.

The integrated development route binds an explicit owner using
`exx_bind_system(system, owner, status)`. The native Hamiltonian, XC/total energy,
Taylor lifecycle, PT-CN freeze/history and restart metadata resolve that owner.
Each owner keeps its own optimized fast-PR kernel, ACE, trial buffers and histories;
there is no process-wide active-owner switch or generic-kernel substitution.
`SALMONPY_HSE_CONTEXT=1` opts the full `main_tddft` driver into this route for
explicit regression runs. Existing binaries and the Python TDDFT API are unchanged.

Scope: serial CPU, periodic scalar/non-spin-polarized HSE, full real-space grid,
full-support exchange and common immutable process settings. The native refresh
and action reject spatial/orbital/k MPI, Wannier, symmetry and spin paths for
bound owners. States can own different cell spacings with the same grid/k counts,
but changing global XC, mesh, propagator or factor settings between calls is
unsupported. Global configuration is a caller precondition, not a fully isolated
configuration object. FDTD still requires common dt.

Keep `owner` alive and never copy/assign it or a bound `system`. Calls must be
serialized. After all calls, use `exx_unbind_system(system)` then
`rt_state_clear(owner)`; clearing destroys the owned FFT plans. A second binding
to an already enabled owner is rejected; clear before reopening. These are
internal Fortran development interfaces, not a stable C/Python multi-instance API.
Native numerical/allocation failures may still stop the process. `rt_state_bytes`
returns status 2 and bytes -1 for resource-enabled owners because opaque FFT-plan
allocation sizes are unavailable; do not interpret this as zero memory usage.

```sh
python tools/test_hse_explicit_state.py HSE_BUILD_RECORD NEW_OUTPUT --native-resources --threads 2
python tools/test_hse_rt_state_runtime.py HSE_BUILD_RECORD REFERENCE_RECORD FIXTURES NEW_OUTPUT --native-context --threads 2
```

The first probe compares interleaved ACE and direct actions for different cells
with separate-process references and checks close/reopen isolation. The second
runs the actual TDDFT driver through 32 synthetic k-mesh/dt/propagator/factor
cases and compares current, energy and prediction/solver traces. Together they
check native HSE ownership plus full-driver numerical regression; they do not
prove interleaving of complete independent TDDFT solvers or physical convergence.
The standalone `test_hse_rt_state.py` applies only to the parent RT-only profile;
use the full-build probes above for the integrated profile.

## Experimental independent HSE resource primitive

`salmonpy/native_adapters/hse_resource_context.f90` is SALMONpy-owned source,
shipped in this repository and wheel. It links to the preceding pinned fast-PR
kernel modules; the full solver sources are also bundled, and compilation remains explicit.
Stage-specific completed files are maintained under runtime/overrides.

`hse_context` owns its FFT kernel and accepted RT/ACE data. `hse_open` fixes the
orthogonal rectangular grid, uniform full k mesh, spacings, screening and hybrid
fraction; `hse_refresh` accepts an occupied orthonormal complex frame, and
`hse_apply` evaluates fraction-scaled ACE action. Arrays use (grid, orbital, k)
order. Source occupations are one per spin, with no spin multiplier; this is an
exchange primitive, not the physical charge-current interface. Source overlaps
must equal the identity within 1e-8 with grid volume weighting. Layout is fixed
strided CPU FFT. It uses the generic kernel action, not the optimized full
adapter dispatch; no performance or full-TDDFT equivalence claim is made.

Call `hse_close` explicitly, including after failed construction; repeated close
is safe. Do not copy/assign contexts, serialize their memory, or invoke them
concurrently. They own opaque FFT plans, unlike copyable RT-only snapshots.
Status 0 means success, 1 invalid input/lifecycle, 2 caught adapter allocation
failure, and 3 a kernel/ACE/numerical failure. Reported refresh failures preserve
accepted ACE/source data, and reported apply failures preserve the output.
Internal FFT work buffers may change. Existing kernel allocation/FFT failures
can still terminate a process; this is not a general failure-atomic public ABI.

```sh
python tools/test_hse_explicit_state.py HSE_BUILD_RECORD NEW_OUTPUT --resource-context --threads 2
```

The probe audits sources and all linked application objects, then tests different
grid sizes, k counts and screening/fraction parameters, interleaved refresh,
independent reference kernel/ACE results, invalid input retention and both close
orders. It is a synthetic ownership regression, not a converged material test.
The primitive is not yet connected to `exx_refresh` or Python TDDFT states.
Wannier, spatial/MPI exchange, symmetry, spin, prediction, total energy/current,
and optimized routing remain outside this stage. FDTD still requires common dt.

## Source ownership and updates

The repository and wheel carry the standard solver and the fast-PR HSE baseline,
including original licenses and provenance. Preparation verifies complete file
inventories, copies sources, and installs stage-specific completed files after
checking both their old and new hashes. It does not use Git or apply patches.
The old UPSTREAM OUTPUT preparation/build arguments remain accepted, but the
bundled pinned source is used and the supplied checkout is never modified.
HSE is supported only through the fast profile; USE_HSE=ON in the standard
profile fails explicitly. Reference-only fast sources are comparison fixtures.

`source_tree_sha256` identifies the bundled tree inventory; `replacements_sha256`
identifies completed stage files. Historical `patch_sha256` values remain solely
for provenance and existing binary receipts, not as evidence of new source
content. NLCC loading additionally rejects unknown source-tree identities when
present. Future numerical source changes require new scientific validation.
Update sources, inventories and their hashes in the same reviewed commit; verify
shared native kernels and rebuild/test affected profiles before distribution.
The optimized native HSE owner is connected as described above; the standalone generic resource primitive remains a separate diagnostic component.


### Shared kernels and field runtime in current fast-PR builds

New explicit-routing fast-PR builds share SBE, FDTD and LD source with the
ordinary profile. Native/solver kernels are byte-checked during preparation.
The same tree builds with HSE ON/OFF; OFF requires no HSE FFTW linkage.

Use NEW_BUILD/build.json for NativeEHField and HSE material groups when the
record declares field_abi protocol1 / salmonpy.eh-centered-current/1. Use the
separate NEW_BUILD/executable-build.json with SalmonExecutable.from_mpi_build
for batch MPI FDTD. Both records bind different binaries by their own SHA256.
Existing embedded field records remain supported; only new capable builds carry
the new declaration. Import never builds or downloads anything.

Field and material groups remain separate processes with host-owned MPI. Full
standard GS/TDDFT source migration and concurrent same-process field/electronic
ownership are unfinished; default preparation commands and old binaries stay valid.

## Experimental full TDDFT contexts with a fixed timestep per state

`salmonpy/native_adapters/tddft_context.f90` provides internal Fortran development
entries `rt_open`, `rt_step`, and `rt_close`. Each context owns the arguments of
the actual native initialization/time-step routines: orbitals, density, Hartree
and XC/local potentials, propagation history, pseudopotential working arrays,
XC handle, HSE owner and output handles. States are initialized independently
from the same restart; initialized contexts are never copied. `rt_step` receives
endpoint A/c and returns charge current density j (the negative of native matter
current) and total energy in atomic units. Every-step energy sampling is required.

The process must first initialize its native configuration and serial runtime;
see `tests/native/test_tddft_contexts.f90`. This is not a stable Python/C ABI.
Scope is periodic unpolarized fast HSE, serial CPU/FFTE, fixed ions, transverse
driving, zero initial A, and common immutable settings except the owned dt. Only the tested
Taylor/ACE, learned Taylor, direct Taylor and PT-CN paths are claimed. MPI,
OpenACC, FFTW Poisson builds, secondary fields, restarts and optional spatial,
projection or decomposed-current output are rejected. Context names use simple
letters/digits/underscore/hyphen and existing output files are not overwritten.

Core settings (grid/k counts, bands/electrons, XC, screening,
temperature, propagator/predictor and factor settings) are captured and changes
are rejected before stepping or initializing another context. This is a guard,
not complete removal of `salmon_global`. Global input initialization still requires a shared configuration; arbitrary global mutation
is unsupported. Different XC/grid/k per context remains future work.
The optional internal `rt_open(..., delta_t=...)` fixes a finite positive timestep
for that context; omission captures the input dt. Initialization, Taylor
coefficients, PT-CN, field evaluation and output times use the owned timestep.
Changing dt during a state lifetime, adaptive stepping and FDTD subcycling are
not supported. Other guarded configuration must remain common.
Do not copy contexts or call them concurrently; use `rt_close` explicitly.
Invalid A or guarded-setting changes are rejected before advancement, but native
allocation/numerical errors may still stop the process. Memory/RSS is not capped.

```sh
python tools/test_tddft_contexts.py HSE_BUILD_RECORD SYNTHETIC_FIXTURES NEW_OUTPUT --reference FAST_REFERENCE_RECORD --threads 2
```

The test advances two distinct A waveforms through the real full TDDFT step,
compares with separate-process references, reverses call order, and closes/reopens
one state while the other continues. It compares current/energy, orbitals,
densities, potentials and HSE histories, and checks returned current sign and
energy against the unchanged native driver using its supported Acos2 pulse.
The native reference rejects arbitrary sampled-field input; it is not patched to
bypass that restriction. Binary traces are synthetic test diagnostics, not
restart/checkpoint files. Local implementation regression does not certify
physical convergence, mixed configurations, MPI or a new Python TDDFT interface.

Use `--different-dt` with the context regression tool to test dt=0.001/0.0015
against separate processes configured with those values. The mixed runs poison
the old global dt and also check reverse order and close/reopen. See
`docs/tddft-state-settings-validation-2026-11-15.md` for local evidence.


## Owned XC workspaces and stencil tables (development)

The explicit fast-HSE profile stores XC grid buffers and PZ work arrays in each
`s_dft_system%xc_work`. `calc_xc(..., workspace=...)` can reuse an explicit PZ
workspace; omission creates call-local storage. Shape/spin changes resize only
that owner. Functional/Libxc handles are not copied or made mutable for scratch.
Modulo tables and the optimized-stencil capability flag are owned by each
`s_stencil` and passed to Hamiltonian/stencil consumers. Context destruction
releases these allocatable buffers automatically.

`tools/test_owned_workspaces.py HSE_BUILD_RECORD NEW_OUTPUT --threads 2` checks
PZ evaluation against the original built-in routines, two shapes, spin/shape
resize, table contents and continued use after the other owner is released.
This is a lower-level workspace test, not heterogeneous full TDDFT. The full
context regression still requires common grid/k/XC: input and other global
initialization settings remain. AVX/SVE tuned kernels,
MPI, GPU and concurrent calls are not validated by the local Apple Silicon test.
See `docs/tddft-workspaces-validation-2026-11-15.md`.


## Owned OpenMP selection (development)

Stencil, current and force OpenMP selection flags are captured in `s_stencil`.
Numerical kernels do not read process-global selection flags. The internal
`rt_open(..., openmp_mode=...)` accepts `auto`, `orbital` or `rgrid`, fixed for
that state and applied to the three kernels. Omission preserves the three
independent input settings and original auto threshold. Invalid values are
rejected before initialization. This selects loop decomposition, not the process
thread count; concurrent host calls remain unsupported.

Add `--different-openmp` to the full-context test for orbital A and real-space B.
Each is compared with its own isolated reference, including reverse order and
close/reopen; mixed runs also poison old input selection strings. It can be
combined with `--different-dt`. Common grid/k/XC restrictions remain. Force
selection is wired but fixed-ion current regressions do not validate ionic forces.
See `docs/tddft-openmp-validation-2026-11-15.md`.


## Explicit construction configuration (development)

`structures:s_rt_run_config` contains `delta_t`, `max_steps`, `stencil_mode`,
`current_mode` and `force_mode`. After native process/input initialization,
`rt_input_config(options)` captures the existing input into a value object.
Set its fields and pass `rt_open(state, name, status, config=options)` to create
an independent state. This copies values; later changes to options do not alter
the state. It does not retain pointers to the configuration.

The configuration requires finite positive dt, capacity 1..huge(default_integer)-1,
and selectors auto/orbital/rgrid. It cannot be combined with the legacy delta_t
or openmp_mode keywords. Invalid/mixed values are rejected before initialization.
The bound prevents endpoint integer overflow, not memory exhaustion. Allocation
can still fail; this is not a memory estimator or a hard resource budget.

RT arrays and last-step handling use the state's capacity; over-capacity steps
are rejected before advancement. `--different-capacity` tests A=12/B=6 steps with
independent dt and OpenMP options, separate references, reverse/close/reopen,
source-options capacity mutation, and global-nt poisoning. The other state's
calculation continues after capacity exhaustion. Grid/k/XC and remaining input
settings still require common immutable values. No general input-object API or
public Python backend is introduced. See `docs/tddft-config-validation-2026-11-15.md`.


## Explicit grid bounds initialization (development)

`structures:s_rgrid_config` holds periodicity (0 or 3), requested point counts,
spacing and consistency-check lengths. `grid_config_from_input(lengths, config)`
captures native input; for Maxwell it preserves the existing al_em/dl_em rule.
`initialization_sub:init_grid_whole_config(lengths, spacing, config, grid, status)`
then reads only its explicit geometry/configuration arguments. The old
`init_grid_whole` remains a compatibility wrapper using that route.

This initializes bounds, halo ranges and index arrays only. Coordinates and MPI
partition maps are left untouched, as in the legacy routine, and must be handled
by the subsequent initialization stages; this is not a complete standalone grid
or solver constructor. Invalid geometry/configuration returns status=1 before
changing the destination. Allocation failure returns status=2; successful trial
arrays are transferred only after all checks. Integer bounds are not a RAM budget.
The legacy wrapper retains process-terminating behavior on failure.

`tools/test_grid_config.py HSE_BUILD_RECORD NEW_OUTPUT` covers even/odd periodic
and isolated bounds, point/spacing inputs, separate shapes, poisoned old globals,
Maxwell consistency and invalid-call retention. This lower-level isolated-grid
test does not add isolated TDDFT support. Cell/k/pseudopotential/XC/HSE initialization
still has shared-input dependencies, so full TDDFT retains common grid/k/XC.
See `docs/grid-config-validation-2026-11-15.md`.


## Explicit cell resolution (development)

`initialization_sub:resolve_cell(config, geometry, status)` uses only the
value-only `s_cell_config` (lengths, vectors, spacing, points) and returns
`s_cell_geometry` (primitive vectors, lengths, spacing, orthogonality).
Orthogonal precedence is explicit lengths, diagonal vectors, then points times
spacing. Tilted cells use their vectors and original norm calculation. Explicit
spacing takes precedence over lengths divided by points. Invalid/nonfinite,
nonpositive or zero-volume geometry returns nonzero without changing the result;
a scaled determinant detects singular cells without forming huge volume products.
This does not certify conditioning of nearly singular cells.

Native DFT/RT initialization calls this resolver; its compatibility layer still
normalizes global al for downstream consumers. This is a local cell resolver,
not completion of global cell/k/pseudo/XC removal or heterogeneous TDDFT support.
`tools/test_cell_config.py HSE_BUILD_RECORD NEW_OUTPUT` checks precedence,
tilted geometry, independence and failed-call retention. See
`docs/cell-config-validation-2026-11-15.md`.


### Owned k-point, pseudopotential and XC inputs

The fast-HSE development profile separates internal Fortran inputs into value types:

- `s_kpoint_config` → `lattice:init_kvector_config(system, unfold, config)`: deep-copies mesh, shift, k-file, unfolding and symmetry inputs. The effective mesh is `system%kmesh`; HSE consumes that mesh and owned symmetry operations.
- `s_pp_input_config` → `salmon_pp:read_pslfile(system, pp, config)`: copies species, files, local channels and masking settings into `pp%input_config`. Detected formats and parser scratch flags also belong to `pp`.
- `s_xc_config` → `salmon_xc:init_xc_config(functional, config, status)`: retains functional identity, spin, coefficients, HSE screening and rVV10 settings. Never copy initialized Libxc handles; call `finalize_xc` before reinitializing. `calc_xc(..., pre_scf_active=...)` takes an explicit evaluation phase, defaulting to production.

Legacy entries remain input-capture adapters. Explicit paths never temporarily swap shared input. These are internal initialization/evaluation components, not new public Python or C interfaces. Specialized SOC/+U PP setup, SCF phase/algorithm controls, communicators and file units still have shared dependencies. Full native TDDFT therefore retains its common grid/k/XC guard. These internal interfaces assume validated input and do not guarantee atomic failure for all native errors or I/O failures.

Run `tools/test_kpoint_config.py`, `tools/test_pp_input_config.py` and `tools/test_xc_input_config.py` with `HSE_BUILD_RECORD NEW_OUTPUT`. See `docs/shared-input-validation-2026-11-15.md` for measured scope and limitations.


### Scalar pseudopotential grid placement

The internal `prep_pp_sub:init_ps(..., config=...)` copies `s_pp_grid_config` into `ppg`.
Periodicity, Poisson/FFT selection and diagnostics are fixed at construction; the legacy entry also retains its first captured configuration.
Atom count/species come from `system`, and species count from loaded PP arrays. Updates with changed atom count, species mapping, grid geometry, radial layout or configuration are rejected. Reuse the same loaded radial numerical tables unchanged: guards check configuration/layout, not table contents, and do not hash or copy full tables each step. Position changes retain the existing `dealloc_init_ps` → `init_ps` sequence and are compared with a fresh grid at the same geometry.
Cache, interpolation and projector communication-table helpers use owned array sizes. NLCC grid placement reads `system%kion` and the periodicity captured by the radial loader.

The explicitly independent path is scalar, without +U. Legacy SOC/+U entries remain, but their specialized shared state has not been extracted. MPI resource ownership, global timers and heterogeneous full-TDDFT restrictions remain. Run `tools/test_pp_grid_config.py HSE_BUILD_RECORD NEW_OUTPUT` for the focused probe.


### Owned XC evaluation phase

The internal fast-HSE entry `exx_set_phase(system, pre_scf)` gives that system an explicit evaluation phase.
`exx_pre_scf(system)` resolves it consistently for HSE enablement, thermal SCF helpers, semilocal XC and rVV10.
RT initialization owns production phase regardless of the legacy flag: another calculation entering pre-SCF cannot change an existing or newly initialized RT state's evaluation.
A phase change clears only the explicitly bound HSE owner’s exchange resources and prediction history. Unbound systems may select a phase before binding, but must not reuse legacy shared exchange caches.
After explicitly returning from pre-SCF to production, refresh exchange from the current orbitals with `exx_refresh` before using its action.

GS/DC pre-SCF controllers now retain threshold, consecutive-success count, required steps and first iteration in each system. GS, orbital solves, DC fragment/total energy and hybrid metadata consult the same phase. The local readiness check, collective reduction, count update, diagnostic, transition and potential rebuild retain their original order. Each SCF invocation resets its count and start iteration; it does not reactivate a completed pre-SCF phase.

Systems without explicit phase ownership and no-argument `exx_enabled()` retain the legacy fallback. DC communicator/material swapping, adaptive exchange controls and legacy exchange caches remain shared. This is not full GS/DC isolation or arbitrary concurrent-call support.

Run `tools/test_pre_scf_state.py HSE_BUILD_RECORD NEW_OUTPUT` for independent controller checks, and `tools/test_pre_scf_runtime.py HSE_BUILD_RECORD ORIGINAL_FAST_HSE_RECORD NEW_OUTPUT` for four tiny serial GS comparisons. The latter uses deliberately loose thresholds: it verifies transitions and scalar outputs, not physical convergence. Add `--dc-ranks 2` or `--dc-ranks 4` and use MPI-enabled build records to compare two DC fragments, each on one or two spatial ranks (`--mpiexec PATH` selects the launcher). The tool disables restart output for spatial hybrid calculations and limits OpenMP/BLAS to one thread per rank. All eight DC/MPI cases matched the original fast-HSE build locally; this does not establish concurrent independent solver isolation. See `docs/xc-phase-validation-2026-11-15.md` and `docs/pre-scf-validation-2026-11-15.md`.

### SCF exchange-control routing

Internal `exx_set_scf_flags(system, ...)`, `exx_adaptive_is_ready(system)`, `exx_is_frozen(system)` and `exx_support_has_changed(system)` resolve the system-bound HSE resource owner. GS SCF uses these for adaptive readiness, temporary exchange freezing and support-event consumption; `exx_check_localization(system=system)` validates that same owner. Collective reductions and freeze/restore ordering are unchanged. Unbound systems still use the legacy owner. This routes existing controls; it does not construct independent full GS solvers, isolate the remaining configuration globals, or enable concurrent calls.

### GS-owned HSE resources

The GS entry now creates a local HSE resource owner and binds it after system initialization when the existing native-owner predicate accepts the configuration. The initial scope is ordinary `theory='dft'`, no geometry optimization, periodic unpolarized non-spinor HSE06, full-support CPU exchange, no Wannier/symmetry reduction, negative temperature setting (integer-occupation route), and one rank in each electronic decomposition. DC, distributed/thermal/localized exchange and other configurations retain the legacy owner. This does not broaden solver physics support.

A selected owner is reported as `EXX_GS_OWNER explicit`; it survives through result/checkpoint output, then is unbound and cleared. `tools/test_pre_scf_runtime.py ... --require-gs-owner` additionally requires that branch in the current build. Full GS state and configuration are not yet independently encapsulated; no concurrent GS API is provided.

The remaining same-input GS lifecycle check is provided by `tools/test_gs_lifecycle.py HSE_BUILD_RECORD ORIGINAL_FAST_HSE_RECORD NEW_OUTPUT`. It explicitly compiles a serial probe, runs two GS solves in one process, and compares both with the original executable. Four pre-SCF configurations match scalar results and every generated restart payload byte, including wavefunctions and density-mixing history. Payloads are deleted after comparison; hashes and logs are retained. This verifies sequential reuse for the tiny fixed-input fixture only, not different-input reuse, concurrent GS contexts, or checkpoint import into TDDFT. See `docs/gs-lifecycle-validation-2026-11-15.md`.

### Retained GS storage (internal Fortran)

`salmonpy_gs_context` provides `s_gs_context`, `gs_run(ctx, status)` and `gs_close(ctx)`. The context retains the former GS work objects: orbitals, densities/potentials, occupations/forces, mixing and pseudopotential data, XC handles and the eligible HSE owner. **Declare actual contexts with `TARGET`, keep their address stable, and never copy a live context.** Close explicitly before discarding one. The external `main_dft` remains a run/close wrapper.

`gs_run` rejects an already live context with status 1 before mutation; successful completion returns 0. `gs_close` is idempotent. Native errors can still terminate the process. Input, parallel and timer globals remain shared; this is not a stable public ABI, an explicit-input GS API, or concurrency support. Only the documented same-input scalar-HSE case is validated for multiple retained contexts. The lifecycle probe compares retained A before/after B's run and close, including electronic arrays, occupations/forces, mixing histories and selected HSE caches; it is not an exhaustive memory comparison.

### Explicit GS SCF controls

`gs_input_config(value)` captures the legacy maximum SCF iterations, convergence metric, convergence threshold, pre-SCF threshold and consecutive-step requirement into `s_gs_run_config`. Pass it as `gs_run(ctx, status, config=value)`; the context validates and copies it before initialization. Omitting it captures the same inputs once for compatibility. `max_iterations` must be positive and below the integer maximum, both thresholds finite and nonnegative, and `pre_scf_steps` positive. Invalid configuration returns status 1 before activating the context.

SCF loop bounds and convergence checks use the owned values, including adaptive-exchange readiness thresholds. Mutating the caller's value or the corresponding legacy globals does not change the selected settings. Occupation-setting inputs, geometry/material inputs and other global controls remain shared. No general different-input/concurrent GS claim follows.

An explicitly positive pre-SCF threshold must also satisfy the original fresh-static-hybrid prerequisites: no restart/MD/optimization, positive temperature for DC, no checkpoint interval/shutdown request, and density mixing/convergence. These are checked again at the context boundary. The legacy executable wrapper reports invalid run configuration as a nonzero failure instead of silently skipping the calculation.

For Python support collecting local files, see [PREPARATION.en.md](PREPARATION.en.md). It does not resolve target compatibility/dependencies or install packages.

The owned `convergence_metric` accepts rho_dne, norm_rho, norm_rho_dng, norm_pot
and norm_pot_dng. SCF arithmetic/normalization formulas and reduction order are
unchanged; selection and labels use the owned value. Supply thresholds in the
internal atomic-unit convention for the chosen metric. Mixing, normalization's
electron-count inputs and general material initialization remain shared.

### Per-state GS mixing settings

Internal `s_gs_run_config%mixing` (`s_mixing_config`) owns `method`, `mixrate`,
`alpha_mb`, `beta_p`, `nmemory_p`, `auto_mixing`, `update_ratio` and
`initial_no_mix_steps`. Capture with `gs_input_config`, optionally edit the value,
and pass it to `gs_run`. It is copied before initialization. Adaptive coefficients,
previous residual and density/potential histories remain separate mutable data in
each state's `s_mixing`. Halving and pre-SCF resets use that state's initial values.

Methods: simple/broyden/pulay/simple_potential/simple_dm; auto_mixing: y/n.
Coefficients and update_ratio must be finite, Pulay memory 1–21, initial no-mixing
steps nonnegative and below the integer maximum. simple_dm still requires
spin-orbit. Storage capacity remains 21; formulas and execution order are unchanged.
This is an internal Fortran contract, not a new Python API or concurrent GS support.
Legacy initialization captures inputs at initialization; the separate DFT-MD entry
still uses a global auto-mixing switch.

Use `tools/test_gs_lifecycle.py ... --mixing pulay --auto-mixing --distinct-mixing`
to compare two states with different mixing settings for the same material against
the original executable under each setting. See
`docs/gs-mixing-config-validation-2026-11-16.md` for evidence and untested scope.

### Electron-count normalization of GS convergence

`s_gs_run_config%normalization` owns `nelec` and `nelec_spin(2)`, captured by
`gs_input_config` and copied by `gs_run`. The original selection is retained:
unpolarized uses nelec; polarized uses the sum of nelec_spin when positive,
otherwise nelec. DC still overrides this with its state-owned `dc%elec_num_tot`.
For non-DC rho_dne, a nonpositive divisor returns status 1 before initialization.
Other metrics and their grid normalization are unchanged.

These are convergence-normalization inputs, **not material electron count,
charge or occupation setters**. Normally retain the captured values. This does
not finish material-input isolation or establish tested polarized/DC operation.
`tools/test_gs_lifecycle.py ... --different-normalization` diagnostically doubles
the second state's divisor: its density residual must halve while electronic
outputs remain unchanged. Evidence and limits are recorded in
`docs/gs-normalization-config-validation-2026-11-16.md`.

### State-owned occupation inputs

Internal `s_dft_system%occupation_input` captures nelec, nelec_spin, nstate, spin, spinorbit, temperature and redistribution_delay at initialization. Initial filling, ordinary ne2mu, the ordinary GS redistribution gate and occupied-band index for gap reporting use that value. Uninitialized systems retain the legacy-input fallback. Temperature uses internal atomic units.

Do not modify a live state's record. This is not a new occupation override for gs_run. The separate DC solver, parallel setup and other HSE settings still reference shared inputs; arbitrary different-count/temperature concurrent GS is not established.

Tests: `tools/test_occupation_inputs.py BUILD_RECORD NEW_OUTPUT` and `tools/test_gs_occupation_runtime.py BUILD_RECORD REFERENCE_RECORD NEW_OUTPUT`. The previous uncommitted files were unavailable after reboot and reconstructed from the conversation. Fresh evidence: `docs/gs-occupation-recovery-validation-2026-11-16.md`.

### HSE state-owned occupation inputs

HSE electron count, temperature and state capacity (nelec/temperature/nstate) now use the same `occupation_input` value. Legacy uninitialized systems keep their fallback. Arithmetic and update order are unchanged.

Independence coverage is serial CPU, periodic scalar unpolarized, integer occupation, full-support fast HSE. Two states with 4/2 electrons are compared in solo/interleaved execution, destruction/recreation and action/exchange energy while legacy inputs are poisoned. Finite-temperature GS is compared with the original solver; independent thermal/Wannier owners are not established. The nstate-dependent Wannier branch is also outside the new independence test.

The outer TDDFT configuration guard remains: other consumers still use shared inputs. This is not a new API permitting arbitrary simultaneous TDDFT electron counts or temperatures. Never mutate a live occupation record.
Run `tools/test_hse_explicit_state.py BUILD_RECORD NEW_OUTPUT --native-resources` (optionally `--threads 2`). Evidence: `docs/hse-occupation-inputs-validation-2026-11-16.md`.

### TDDFT electron-count diagnostics

The `subdip` electron-count check uses the state's `occupation_input.nelec/nelec_spin`. The existing precedence is preserved: for two spin channels, use the positive sum of spin counts; otherwise use nelec. Strict inequalities, output and termination behavior are unchanged. Uninitialized systems retain legacy-input fallback.
This is a coarse anomaly check (strictly between 0.1 and 10 times the specified count), not a precision conservation or physical-convergence guarantee. Shared I/O, ionic inputs and outer TDDFT guards remain.
Test: `tools/test_occupation_inputs.py BUILD_RECORD NEW_OUTPUT --dip`. Scope: `docs/dip-occupation-inputs-validation-2026-11-16.md`.

### Direct fast-PR HSE checkpoint (internal Fortran entry)

`native_adapters/tddft_context.f90` provides `rt_checkpoint(ctx,path,restore,status)`: false writes numerical state, true loads it immediately after `rt_open` from the same original GS/input in a fresh process. Continue at `rt_accepted_step(ctx)+1`; no elapsed time evolution is repeated. Rebuild static exchange/FFT resources in the new process; never serialize raw FFT handles.

The payload covers both orbital buffers, density/local/Hartree potentials and histories, RT history, fields/current/energy reference/accepted step, ACE, learned-factor/PT-CN histories and exchange caches. The caller must authenticate original inputs/GS, runtime/adapter and archive hashes before loading. This same-binary native format is not a generic portable archive. Settings mismatch is rejected; a late restore failure invalidates further stepping.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python tools/test_hse_direct_checkpoint.py FAST_HSE_BUILD_RECORD SYNTHETIC_FIXTURES NEW_OUTPUT
```

The tool explicitly compiles adapters against audited fast-PR objects and compares separate continuous/save/resume processes. Two/three k points, Taylor ACE/learned ACE/full-action within the fast-PR engine, and learned PT-CN, split6/7, match currents,total energy and internal traces exactly after6/5 remaining steps. Validated at OpenMP1/2. Conventional HSE and comparison baselines are not execution engines.

Scope: internal serial CPU Fortran entry, periodic scalar fixed occupation/ions, full support. This does not provide a stable native ABI, GPU, symmetry, orbital/grid-distributed ACE or reduced-q support. Explicit host-MPI k-only HSE and FDTD coupling use the separate entry below. Synthetic numerical regression is separate from physical convergence. Historical native text-file positions are not resumed.

Save at an accepted step with initialized exchange resources. Step-zero saving immediately after initialization is not established.

For outer-MPI/FDTD fast-HSE coupling, explicitly build the dedicated worker and use FastHSEResponseBackend. In this worker mode each outer material rank calls its own serial process. See COUPLED.en.md for build/coupling/direct-restart instructions. This does not add inner HSE MPI.

For internal k-only MPI use `tools/build_fast_hse_mpi.py` and `salmonpy.fast_hse_mpi.NativeFastHSEBackend`. Supported material groups have one or two ranks, no orbital/rgrid split, and one native lifetime per process on the main thread. The host owns MPI; the backend temporarily changes cwd. Same-topology rank-local direct restart preserves numerical HSE histories and rebuilds static FFT resources. See COUPLED.en.md for the four-rank FDTD/SBE/HSE example and failure limits.

The fast source now shares ordinary initialization input storage and RT output
management too. The duplicated adapter output source is removed. Existing MPI
and serial builders reuse the compiled solver helper, or compile the same shared
source for ABI-compatible older prepared builds. HSE-specific settings retain their own complete
snapshot; ordinary retained GS/embedded TDDFT migration is still unfinished.

Fast HSE model registration now uses the same model-owned parser workspace
contract as the standard profile: scratch connections, namelist diagnostics,
coordinate-source flags and common input/default snapshots belong to each read.
The full fast-HSE defaults, conversion rules and checks are retained; registered
models still capture all HSE-specific settings in hse_model_inputs. The common
parser values record alone is not a complete HSE model snapshot. Legacy batch
reads keep their named temporary files for MD compatibility. This is parser
migration, not completion of retained GS/ordinary embedded TDDFT integration.

The fast native GS context now has an internal import-only path:
`gs_run(..., import_only=.true.)` reads an existing compatible GS wavefunction
restart through the original initialization, without SCF, PBE warmup or GS
checkpoint output. States retain independent XC/HSE/MPI resources; closing one
must preserve the others. `gs_export_seed` returns detached orbitals, occupations,
common input values and ownership/grid metadata using the shared gs_seed module.
This seed is not a complete HSE model or portable checkpoint. Python
NativeGSBackend.load now accepts a fast build explicitly declaring gs_abi/1.
That profile returns density/occupations/orbitals and material IDs; energies,
residuals, tolerance and convergence remain None. Standard-runtime SCF behavior
is unchanged. Fast session.run, eigenvalues/forces and archive export remain unconnected.
New builds additionally declare gs_abi.memory_rt=1 after probing the compiled
memory GS protocol: retained GS→fast-PR HSE model creation is supported without
GS restart files/SCF/replay (NATIVE_GS.en.md). Ordinary XC is supported by the HSE-disabled fast profile described below. The imported
state does not certify convergence, energy or forces. MD/DC/optimization import
and real-material convergence validation are outside this measured path.

The fast source also builds an ordinary-XC host library via build_embedded_tddft.py with USE_HSE=OFF. It declares source_profile=fast-pr, ordinary_rt_abi=1 and import-only gs_abi/1. Shared embedded RT/worker/parallel-bridge sources preserve normal APIs; fast numerical states already own stencil/XC scratch. PZ/TB-mBJ/Libxc-PZ/PBE memory GS handoff, independent states and direct RT restart are supported as described in NATIVE_GS.en.md. Computed GS and export migration remain pending.
