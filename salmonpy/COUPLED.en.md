# MPI TDDFT coupled to Multiscale1D or native E/H FDTD (experimental)

[日本語](COUPLED.md) · Topic: `coupled` · Paired revision: `2026-11-19.1`

Split one Python MPI job into a field group and material groups. The legacy/default
selection uses one material and one macro point per material group; explicitly
shared or registry routes described below can own multiple points/models per group.
Every member participates in that group's native electronic calculation. k-point,
orbital and real-space decomposition stays inside the group. There is no subprocess
or nested MPI launch. Native FDTD retains many local LD states.

The entries are `NativeTDDFTMaterialProvider`, `NativeTDDFTBackend`, `NativeEHField`
and `CoupledStepper`; the latter two live in `salmonpy.coupled`. Existing serial-worker
TDDFT and batch FDTD remain available with their own restrictions.

## Field runtime records

NativeEHField accepts the existing embedded-tddft/1 host-library record and new
fast-hse-mpi/1 records declaring field_abi protocol1 with the same E/H timing
scheme. A new build_fast_hse_mpi build can therefore supply both FIELD_BUILD and
HSE_BUILD in the mixed_hse_models example; rank groups stay separate processes.
Its executable-build.json is the separate batch MPI FDTD executable record.
This does not add standard TDDFT/GS support to a fast library merely by changing
the record schema, and import does not compile or acquire a runtime.

## Choosing dimension and response model

Three-dimensional fields and macro-point layouts are possible, but TDDFT keeps an
electronic state at every macro point. **Start with one-dimensional propagation or
layered structures for practical TDDFT work.** This refers to the electromagnetic/
macro layout, not reducing the microscopic DFT/TDDFT electronic problem to 1D.

**For three-dimensional calculations with electronic feedback, recommend SBE or a
user-written current solver.** Three-dimensional TDDFT is allowed, but choose it only
after assessing the required physics, independent state count, memory and runtime.
SBE/custom models must still satisfy their physical validity range, independent
point-state ownership and A/c-to-charge-current-density contract. The ordinary mixed
sample connects SBE to native E/H; dimension selection remains explicit.
Standalone three-dimensional LD FDTD is a separate workflow.

The current 40³-grid/two-TDDFT-point sample is a connection regression, not a general
3D TDDFT recommendation or a transversely uniform 1D multilayer implementation.
Simply reducing transverse grid counts to one does not implement a 1D solver.

## Preparation and resources

First establish intent, total memory allowance and CPU concurrency. The standard
sample uses two field ranks plus two two-rank material groups: six ranks and two
electronic states. Start with one OpenMP/BLAS thread per rank. The mixed sample offers uncertain memory estimates and measurements; excess only warns. Measure a small GS before scaling.

Explicitly build a **shared library against the host MPI**. Portable serial runtimes
and the communication-only bridge cannot substitute for it. Nothing builds on import.
The current build route targets macOS/Linux and GNU-compatible compilers.

```sh
python tools/prepare_salmon_runtime.py UPSTREAM_CHECKOUT NEW_SOURCE
python tools/build_embedded_tddft.py NEW_SOURCE NEW_BUILD --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=mpifort \
  --cmake-arg=-DCMAKE_C_COMPILER=mpicc \
  --cmake-arg=-DBLA_VENDOR=OpenBLAS \
  --cmake-arg=-DOPENMP_FLAGS=-fopenmp \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --cmake-arg='-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -ffree-line-length-none -fallow-argument-mismatch'
```

Resolve compiler/BLAS locations for the host. The C wrapper must accept the selected
OpenMP flags too. The tool enables MPI/PIC, audits sources and objects for forbidden
eigensolvers, and records library SHA256, source hashes and commands in `build.json`.
Tools also ship under `salmonpy/_runtime_source/tools`; locate them with
`importlib.resources.files('salmonpy').joinpath('_runtime_source')`.

mpi4py and the library must use the same MPI. This host-specific library is not a
portable wheel runtime. Request a Libxc-enabled build explicitly if needed; current
inner-decomposition numerical validation covers PZ and TB-mBJ.

## Sample

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 6 python -m salmonpy.examples.mpi_coupled_tddft \
  NEW_BUILD/build.json GS_BUNDLE RT_INPUT NEW_OUTPUT
```

The self-contained `GS_BUNDLE` contains saved GS records, relative pseudopotentials
and the RT input named by `RT_INPUT`. Match GS functional/grid/settings. The loader
does not run DFT or convert restart formats. `NEW_OUTPUT` must be a new location on
a filesystem shared by all members of each group.

The sample uses a 40³ field grid and two adjacent electronic macro cells with an
artificial pulse. It is a connection example, not converged film optics. Using
`--field-ranks 1 --material-ranks 2` requires five total ranks. Inner decomposition
is `(k, orbital, x, y, z)` and its product must equal material-group size; choose
counts compatible with the actual k points, bands and real-space grid. The sample
splits orbitals.

Providers and response backends in `MaterialBinding` remain replaceable user objects.
Loading/preparing material data is separate from A-to-J physics. Return physical
**charge-current density**, avoiding double conversion of electronic-current sign.
Group collectives use the supplied `solver_comm`.

## One-rank Multiscale1D with distributed TDDFT

`CoupledStepper.from_multiscale(context, bindings, field=None)` connects a Python
`Multiscale1D` state to the same material groups. Construct the field state with
`model.new_state(dt=dt, steps=steps)` on the single field rank; pass `None` on
material ranks. All ranks call the factory and then `step()`. Material cells
1 through N map to `context.plan.point_ids` in that order. Two independent
material groups may prepare the same GS; each owns its own evolving state.
One field rank plus two two-rank material groups requires five MPI ranks.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 5 python -m salmonpy.examples.mpi_multiscale_tddft \
  NEW_BUILD/build.json GS_BUNDLE RT_INPUT NEW_OUTPUT --material-ranks 2 --steps 400
```

Use `--material-ranks 1` with three total ranks. The sample uses two material
cells, spacing 5 a.u., vacuum (300, 50), dt=0.02 a.u. and a z-polarized artificial
pulse. Its orbital decomposition must fit the supplied GS. `response.json`
contains endpoint A/j histories and time declarations, `field.npz` contains the
whole-grid A history including t=0, and `communication.json` contains measurements.
Outputs are not restart files. The sample does not build or prepare a GS.

The factory checks material-site A(0) before creating native states. A nonzero
value is rejected with its point/value: increase the left vacuum cells or delay
the incident pulse. Incident A in vacuum may be nonzero. Do not silently zero
material A(0) or substitute an initially illuminated GS. This first connection
supports only zero initial material A and one native lifetime per process.

The update remains previous j → Multiscale1D field step → new A → material step.
Initial feedback is zero. Both clocks are at n·dt (offset 0); the E/H half-step
adapter does not apply. `step()` returns point-ordered `A_over_c` and `J_charge`
on the field rank. Keep the original field state to read its `grid_A` and
`grid_E`; the latter retains its centered-difference time (n−1)·dt.
Each returned frame declares `E_time`: Multiscale1D's centered E is at
(n−1)·dt while A is at n·dt; native E/H reports E at n·dt. Do not align E and A
by array index alone. `coupled.response.initial_current` is a point-ID dictionary
of physical t=0 currents; a custom state that does not expose one gives `None`.
This observable is separate from zero initial field feedback. The sample records
it as `initial_current` in `response.json`.

Close `coupled.response`, then `coupled.field` where present, then the Python
context. A failure after evolution invalidates the coupling; this is not a
cross-group atomic transaction or a restart checkpoint.

Call `coupled.response.communication_report()` collectively for per-rank
communication measurements. A distribution and current collection are measured
separately from material execution; these elapsed times can include waiting for
other ranks and must not be interpreted as pure network latency. A timing failure
is a measurement failure, not evidence of a failed calculation. Report unavailable
measurements explicitly. Single-host timings do not establish multi-node scaling.

Acceptance compares full-grid A and each material point's j against two-point
SALMON Weyl multiscale output, alongside serial-worker and rank-count parity.
A synthetic GS checks the interface/update order, not physical convergence.

## Native E/H temporal and spatial contract

The new scheme is `salmonpy.eh-centered-current/1`. It is a different discretization
from the Weyl multiscale path; bitwise equivalence between those schemes is not
claimed. Native E/H ordering and incident-source evaluation times are retained.

| Quantity | Field clock |
|---|---|
| E after update n; native saved E/H (H is time-averaged) | n·dt |
| New A and J; raw updated H | (n+1/2)·dt |
| Previous J consumed by E update | (n−1/2)·dt |

Electronic local time zero corresponds to field time dt/2. Initial E, H, A and
feedback J are zero. A nonzero computed initial GS current does not replace that
initial feedback. The sequence is previous J → native E/H update → A update →
TDDFT J update. `A_new = A_old − dt * S(E_n)`; crystal drive is always A/c.

This table is declared in `salmonpy.time/1` form (`NativeEHField.timing`: j in at
offset ½ with lag 1 and t=0 unused, A out at offset ½). `CoupledStepper` gathers
every point's matter declaration over MPI, places it on the electronic clock
(field clock + ½, `ELECTRONIC_CLOCK_OFFSET`) and rejects mismatches before the
first step; `stepper.timing` and `response.json` record both. See
[FIELDS.en.md](FIELDS.en.md).

S matches native observation collocation: Ex averages edges `(i,j,k)` and
`(i−1,j,k)`; Ey/Ez use the corresponding previous y/z edge. Deposit J using S transpose,
half to each corresponding edge. Each macro point represents one uniform grid-cell
volume, so no extra volume multiplier applies. Each field rank adds only to edges it owns.

Coupled points stay at least ten cells from boundaries and clear of the incident
source plane/support and nonvacuum constitutive edges. The surrounding 3³ background
cells must be vacuum. LD may occupy other regions.

## Native E/H state, outputs and limitations

- One native state lifetime per process. No reinitialization, multiple states in
  one process or thread sharing. Two groups may evolve independent copies of the
  same GS. Start a fresh MPI job for another trial.
- Input, initial state, dt and decomposition are fixed at creation; only A or
  previous J varies per step.
- Calls temporarily switch the process working directory for legacy relative I/O
  and restore it afterward. Main thread only; do not perform concurrent relative
  file I/O in the same process.
- Advancing a field/material independently after coupling is rejected as lost
  synchronization. Coupled-step errors still invalidate the session. Standalone TDDFT numerical failures have the scoped recovery contract below.
- Collectively close response, field, then Python context. Fortran does not finalize
  MPI. Native hangs/failures require launcher timeouts; the sample aborts the job on error.
- Same-topology protected restart is documented below. Nonzero initial fields,
  impulse sources, rotations, unequal macro volumes, subcycling and TTM remain unsupported.

`response.json` stores point-ordered A/c, sampled E, previous/new charge currents
and both clocks. `salmonpy-coupled-field.json` records mapping, MPI decomposition,
library and completion status. Material work copies/native outputs remain available.
Select large spatial outputs and sampling intervals in the RT input. These JSON
files are not complete restart checkpoints.

Local regressions check original vacuum/LD FDTD parity, current injection at MPI
boundaries in all three directions, one/two field-rank parity and TDDFT A-history
replay. Physical grid/time/GS convergence and production scaling remain separate work.

## Different materials and time steps

All material states connected to the same current FDTD loop use the field's common dt, even when their material, grid or XC differ. Independent states having separate dt values does not imply multirate coupling support. The coupler checks every point's timing declaration before stepping and rejects mismatched dt. An integer ratio, such as a material step half the field step, is also rejected without a separately validated adapter.

Preserve the native Yee half-step clock and lag1 current feedback. Equal dt does not permit a mismatched half-step phase. Subcycling, intermediate A interpolation and time-averaged current feedback are not implemented or applied automatically. Different dt requires independent runs or a future explicitly validated temporal adapter.

## Toward mixed LD, SBE and TDDFT

`mpi_mixed_multiscale` couples LD and a synthetic SBE point to one 1D field:
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_mixed_multiscale NEW_OUTPUT
```
One field rank and two material ranks check nonzero responses and exact agreement with an independent serial schedule. Output includes point IDs10/20, material keys, A/c, charge current j and communication timing. Endpoint LD derives E from the A difference; it is not the native E/H semi-implicit LD scheme.

With four ranks, `--runtime RUNTIME --bundle GS_BUNDLE --input rt.inp` prepares a third managed serial TDDFT worker point. Four-rank,100-step three-way execution using the Si starter GS (PZ,16³ grid,2×2×2 k mesh) matches an independent serial calculation exactly. Supply a compatible GS/runtime explicitly. This validates the managed serial worker, not internally MPI-parallel TDDFT or converged optical properties. Evidence: `docs/mixed-three-material-validation-2026-11-16.md`.

The final topology keeps many LD states local to FDTD ranks while electronic groups own SBE/TDDFT states. Do not copy the demonstration's one-LD-cell-per-rank layout into large production jobs. This is not real-material convergence or validated 3D mixed coupling.

## Field-local native LD with SBE and internally parallel TDDFT

`mpi_mixed_native` uses one field rank holding two LD regions, one SBE rank and two ranks for one TDDFT point. Four outer ranks exercise orbital decomposition inside the TDDFT group.
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT
```
Use three outer ranks with `--tddft-ranks 1`. `--without-ld` replaces LD regions with vacuum. BUILD_JSON is an explicit `tools/build_embedded_tddft.py` build for the host MPI; SERIAL_RUNTIME provides an independent response replay. Nothing builds or downloads implicitly.

LD retains native E/H coefficients, auxiliary variables and update order. Electronic points retain A/c drive and charge current j. Point10 is SBE,20 TDDFT; constitutive IDs1/2 belong to the field material grid. Electronic neighborhoods do not overlap LD. Field half-step clocks are distinct from electronic local time.

This tiny40³/2-electronic-point 3D fixture is a connection regression; keep recommending 1D propagation for practical TDDFT. It does not establish large-scale LD performance, multiple distributed field ranks, whole-session restart or physical convergence. Native failures are not recoverable; the sample aborts its MPI job.

Validated:128 steps,1/2 internal TDDFT rank field/current comparison, serial response replay, vacuum-LD contrast, original LD-only parity and installed-wheel execution. Scope/errors: `docs/native-mixed-material-validation-2026-11-16.md`.

### Parallelize the field group too

Two field ranks along y, one SBE rank and one TDDFT rank:
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT --field-grid 1 2 1 --tddft-ranks 1 --compare REFERENCE_OUTPUT
```
REFERENCE_OUTPUT is the same sample with one field and one TDDFT rank (three outer ranks). Omit `--compare` to run independently. Comparison requires matching steps,dt,material configuration and point order, then checks A/c,E,j histories. Required outer ranks: `product(field-grid)+1+tddft-ranks`.

Two-way x/y/z field decompositions each matched the one-field-rank result exactly over128 steps. The y/z layouts split both LD materials across both field blocks. Simultaneous two-field/two-TDDFT-group ranks require five ranks and were not tested within this four-core campaign. This is not a large-scale performance or memory-cap guarantee. Evidence: `docs/distributed-field-mixed-validation-2026-11-16.md`.

## Assemble ordered layers before launching MPI

`salmonpy.layers.assemble_native_layers` accepts a `salmonpy.native-layers/1` layout and an explicit mapping from material IDs to LD models or `MaterialBinding` objects. Use existing configured material loaders/bindings for physical parameters and initial states; the layout references their IDs. Custom response backends use the same `MaterialBinding` contract.

```python
import json
from salmonpy.layers import assemble_native_layers

layout = json.loads(layout_path.read_text())
assembly = assemble_native_layers(
    layout, materials, grid=grid, dt=dt, steps=steps, source=source,
    field_grid=(1, 1, 1),
    ranks_per_point={"sbe-point": 1, "si-tddft": 2}, rank_budget=4,
)
report = assembly.record()
print(report["required_ranks"], report["layers"])
# Pass assembly.parallel_plan to ParallelContext.create, assembly.calculation
# and assembly.points to NativeEHField, and assembly.bindings to ResponseSession.
```

The variables above are explicitly prepared models, Grid, PlaneWave and calculation settings, not automatically downloaded assets. Planning starts no solver or MPI process. Launch exactly the reported rank count with the site's MPI launcher.

`examples/mixed_layers.json` provides an ordered SBE/Si/vacuum/Drude/Lorentz fixture. Lengths are explicitly Bohr. Each layer has a unique `id`, a `material_id` (null for vacuum), `thickness` and two `transverse_cells` counts. `axis` must match the source; `origin_cell` gives the first axial cell and transverse anchor. Thickness rounds upward by the axial grid spacing. The report preserves requested/actual thickness, exclusive cell bounds, material identity and ranks. Even transverse counts are anchored at `anchor-count//2`.

Every electronic grid cell gets its own state and MPI group; LD cells stay on field ranks. These are finite rectangular regions, not automatically homogenized infinite films. Increasing transverse area or thickness can greatly increase electronic ranks. The planner warns on an exceeded rank budget and preserves every point and the required rank count. Execution still requires that many MPI ranks. It rejects unknown IDs, incompatible backends, absorbing-boundary overlap and electronic points in an LD constitutive neighborhood. It inserts no implicit gap. The rank budget is not a memory limit (`memory_bytes` is unknown). Initial-state preparation, MPI startup and physical convergence remain explicit.

Run the installed fixture by adding `--layers PATH_TO_mixed_layers.json --rank-budget 4` to the four-rank `mpi_mixed_native` command above. The sample recognizes `sbe-point`, `si-tddft`, `drude` and `lorentz`; the public planner accepts arbitrary explicitly bound IDs. An SBE/LD-only layout needs no TDDFT runtime or GS; the sample's positional runtime/bundle arguments may then be unused placeholders. Pure LD uses FDTD directly.

### Independent points sharing one material

Read the packaged layout and increase the Si thickness from2 to4 Bohr to create two independent TDDFT points from the same initial state:

```python
import json
from pathlib import Path

layout = json.loads(Path("mixed_layers.json").read_text())
for layer in layout["layers"]:
    if layer["material_id"] == "si-tddft":
        layer["thickness"] = 4
Path("two_si_points.json").write_text(json.dumps(layout, indent=2))
```

Run `mpi_mixed_native` with four ranks and `--layers two_si_points.json --tddft-ranks 1 --rank-budget 4` (field1 + SBE1 + two TDDFT groups of one rank). Material ID is shared; point IDs, MPI groups and evolving states are distinct. Increasing thickness also shifts subsequent layers; inspect reported cell bounds.

To repeat SBE instead, change the condition above to `sbe-point`. Both layouts were checked over128 steps by independently replaying each point's distinct A/c history and comparing currents. This does not test multiple native TDDFT states within one process.

### Replay restart for the coupled sample

`mpi_mixed_native --save-at 64` stops at step64 and writes input fingerprints and A/E/current histories to `OUTPUT/checkpoint/`. Use the same total steps, layout, library, GS and MPI topology with a new output directory and `--resume PREVIOUS_OUTPUT/checkpoint`. Add `--compare CONTINUOUS_OUTPUT` to compare complete histories.

This reconstructs field, LD and electronic histories by rerunning their original schedule from initial conditions. Recorded A/E/current must match exactly before continuation. It requires the original assets and repeats the elapsed computation; it is not a direct internal-state dump. Scope is the supplied sample in the same environment/topology. General coupled-state archives remain unimplemented; see the direct sample restart below.

### Direct internal-state restart

Use `mpi_mixed_native --direct-save-at 64` to stop at64 and save rank-local state in `OUTPUT/direct-checkpoint/`. Start a new output/job with the same inputs, library, GS and MPI topology and `--direct-resume PREVIOUS_OUTPUT/direct-checkpoint`. Keep the original total128 steps. No propagation to the saved time is repeated. `--compare CONTINUOUS_OUTPUT` compares histories; `result.json` reports `restored_step` and `native_steps_this_job`.

Saved state includes native E/H/PML, LD P/J/Q and old fields, electronic orbital buffers/density/potentials/Hartree/RT histories/A/current, SBE density matrix, and Python integrated A/lagged current/clocks. Initialization rebuilds static coefficients/resources before loading dynamic arrays. A completed file inventory with SHA256 checks is required; changed inputs/library/topology are rejected.

Scope is periodic middlepoint propagation with the same environment/decomposition: PZ, TB-mBJ (automatic/fixed cval), Libxc-PZ/PBE. Libxc requires an explicitly built compatible native library. Generic custom-backend persistence, changed MPI partitions and other XC profiles require further implementation/validation; fast-PR HSE has a separate entry below. Protected candidate restore is now the default and preserves originals on handled staging failures; the explicit unprotected route has no rollback guarantee. Allow storage and memory for full states.

Direct checkpoint format v2 validates XC names, derivative settings and input cval while retaining freshly initialized Libxc handles. Preserve runtime/checkpoint hashes. The sample RUNTIME argument also accepts an explicitly built/audited serial SALMON executable for reference validation; do not use an MPI binary as this child process. Use matching GS and RT functional inputs for each XC.

#### Internally parallel TB-mBJ/Libxc-PBE direct restart

Use matching-XC GS/RT bundles and four outer ranks: field1,SBE1,TDDFT orbital ranks2. Supply `--tddft-ranks 2 --rank-budget 4 --layers LAYOUT_JSON`. Run a continuous reference, a save job with `--direct-save-at 63`, and a fresh output/job with `--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT`. Keep total steps128.

For TB-mBJ (automatic coefficient/fixed cval1.2) and Libxc-PBE, the65 resumed steps produced complete A/c,E,current histories bit-identical to continuous128. This validates direct same-topology two-level MPI restart. Changing decomposition remains unsupported. Field2+SBE1+TDDFT2 needs five ranks and remains outside this four-core campaign.

Checksum verification on direct restart is partitioned across ranks. Every rank must read the same completed manifest and all checks must succeed before state restoration. Payloads must reside in one shared directory visible with identical contents to all ranks; different node-local copies are outside this contract. `result.json` records per-rank file/byte counts in `checkpoint_validation`.

## Fast-PR HSE in outer MPI / FDTD coupling

`FastHSEMaterialProvider` / `FastHSEResponseBackend` connect a dedicated serial fast-PR worker to the existing A/c→charge-current j contract. FDTD keeps its native MPI entry, LD stays on field ranks, and each HSE point's single-rank group owns an independent child process. This is outer MPI coupling, not internally parallel HSE. Host cwd is unchanged; import starts no MPI, build or download.

```python
from salmonpy.fast_hse import FastHSEMaterialProvider, FastHSEResponseBackend
from salmonpy.material_response import MaterialBinding

binding = MaterialBinding(
    "hse-point",
    FastHSEMaterialProvider(bundle_directory, "rt.inp", max_steps=256),
    FastHSEResponseBackend(worker_build_json, existing_work_root, threads=1),
)
# Pass this binding to existing assemble_native_layers / ResponseSession.
```

Explicitly link the worker from an audited fast-PR build into a new external directory. Conventional HSE is not an engine option.

```sh
python tools/build_fast_hse_worker.py FAST_PR_BUILD_JSON NEW_WORKER_DIRECTORY
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_mixed_native FIELD_BUILD_JSON UNUSED_REFERENCE HSE_BUNDLE NEW_OUTPUT --fast-hse-build NEW_WORKER_DIRECTORY/build.json --layers PATH_TO_mixed_hse_layers.json --dt .001 --steps 256 --rank-budget 3
```

The packaged layout contains two LD materials,SBE,HSE. `UNUSED_REFERENCE` is an unused positional argument for this HSE-only example and may be a nonexistent placeholder; independent reference validation uses another identical worker. For two field ranks use `mpiexec -n 4 ... --field-grid 1 2 1 --rank-budget 4`. Each HSE group requires one rank. The child uses that group's CPU while its MPI caller waits. Budget extra CPUs explicitly before increasing OpenMP threads.

Supply an original GS/RT bundle with tddft_response,hse06 and an explicit fast propagator (hse_taylor4/hse_taylor4_full/hse_ptcn). Only selected dt/step capacity are overridden; other physics inputs remain. Use serial internal parallel settings and the fast engine's dense-k/fixed-occupation/full-support restrictions for learned modes. The zero-initial-field sample requires initial HSE current<=1e-12 au and subtracts no background current. The preparation tool generates a synthetic continuation fixture, not a converged material GS.

`--direct-save-at 127` captures field/PML/LD/SBE and worker internal state at one epoch. A new job uses `--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT`. LOAD restores accepted step/current, with saved A/coupling histories checked. Only129 remaining steps are propagated; no elapsed127-step reconstruction. Fresh static initialization and independent serial validation replay remain separate. Same worker,bundle,input,MPI topology and shared archive directory are required.

POSIX macOS/Linux worker uses /dev/stdin; measured here on Apple Silicon. Orbital/grid-parallel HSE, changed-topology restore,GPU,generic backend archives and continuation after native failure remain unsupported.

Worker inputs must use atomic units (unit_system=au), no internal pulses: ae_shape1=none or explicit zero impulse, ae_shape2=none, trans_longi=tr. Double driving by an internal pulse and external FDTD A is rejected. State time is accepted step×dt.

## Fast-PR HSE with internal k-point MPI

`salmonpy.fast_hse_mpi.NativeFastHSEBackend` runs fast-PR HSE collectively on an explicitly supplied material communicator. The host initializes MPI and retains ownership of finalization. Import does not load MPI or compile. Each group owns one point and one native lifetime per process, on the main thread. The native call temporarily changes the process working directory: do not use other threads that access relative paths.

Prepare a new explicit-routing source with `tools/prepare_hse_rt_state.py NEW_SOURCE --explicit-routing`, then explicitly build it with `tools/build_fast_hse_mpi.py PREPARED_SOURCE NEW_OUTPUT --jobs 2 --cmake-arg=...`. Use matching host-MPI Fortran/C compilers, OpenBLAS/LAPACK, FFTW and optional Libxc paths as in the runtime build guide. This build requires MPI/HSE and rejects GPU, Poisson-FFTW, ScaLAPACK and EigenExa profiles. The receipt binds source and adapter hashes. Only the pinned fast-PR profile is an engine option.

```python
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend
from salmonpy.material_response import MaterialBinding

binding = MaterialBinding("hse-point",
    FastHSEMaterialProvider(bundle_directory, "rt.inp", max_steps=256),
    NativeFastHSEBackend(mpi_build_json, existing_shared_work_root))
# Assign two ranks to this point with assemble_native_layers / ParallelPlan.
```

Supported decomposition is k-only, one or two ranks, `nproc_ob=1,nproc_rgrid=1,1,1`. The backend overrides `nproc_k` to the communicator size and overrides dt/step capacity explicitly. Use enough k points for the group. Other driving/fixed-occupation/full-support restrictions above remain. The sample requires a serial worker receipt too, for independent reference validation after the coupled run.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native FIELD_BUILD_JSON UNUSED_REFERENCE HSE_BUNDLE NEW_OUTPUT --fast-hse-build SERIAL_WORKER_BUILD_JSON --fast-hse-mpi-build MPI_BUILD_JSON --hse-ranks 2 --layers PATH_TO_mixed_hse_layers.json --dt .001 --steps 256 --rank-budget 4
```

This assigns field1,SBE1,HSE2; both LD materials stay local to the field rank. Use three separate jobs with the same arguments: continuous, `--direct-save-at 127`, then `--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT`. Rank-local HSE shards include group size/rank and k ranges; all shards are checked before committing numerical state. Fresh initialization recreates static exchange/FFT resources, not elapsed propagation. Restore requires the same binary, bundle, communicator topology and shared archive directory. Early invalid/truncated shards are rejected before numerical commit. Default protected restore preserves original states on handled staging failures. Explicit unprotected restore and unhandled native failures have no retry guarantee.

Measured on Apple Silicon with a synthetic 8³ microscopic grid and 4×1×1 k mesh: Taylor/learned PT-CN, field1+SBE1+HSE2, 256 steps at dt=.001, direct split127. All A/c,E,j history arrays matched bit-for-bit after129 resumed calls. Against independent serial fast-PR response, maximum HSE current differences were4.07e-16 and1.71e-16 au. These are connection/restart regressions, not material convergence or scaling certification. Evidence: `docs/fast-hse-kmpi-validation-2026-11-16.md`.

During layer planning, backends may declare `supported_group_sizes=(1,2)`: a nonempty tuple of unique positive integers. Unsupported rank counts are rejected before MPI or initial-state preparation. Fast HSE declares `(1,2)`. Custom backends without this optional declaration retain existing execution/decomposition checks. This is not certification of numerical/restart capabilities and ResponseSession also validates it and agrees declarations across all groups before preparing materials.

The direct sample writes manifest schema `salmonpy.native-direct/2` with persistence format,save/restore modes,group rank/point ownership and field/electronic k-orbital-real-space decompositions. All ranks compare it with the current candidates before any LOAD, including layout changes with the same total rank count. Static resources initialize first; preflight performs no propagation,snapshot,restore or LOAD. Legacy manifest/1 lacks this metadata and is rejected by the new route; use the original compatible reader for old archives. Native payload formats are unchanged. Custom direct-state/1 callbacks require explicit checkpoint_format,checkpoint_modes (tuple of save/restore) and checkpoint_decomposition (positive-integer tuple whose product equals group size). Declarations are caller contracts, not certification of implementations or atomic rollback. A snapshot method alone does not establish direct restart support; the supplied helper accepts only built-in LD/SBE snapshot routes, not replay-based states.

Capability and MPI compatibility checks use the common `salmonpy.checkpointing` API. See [CHECKPOINTING.en.md](CHECKPOINTING.en.md) for custom backends and serial hosts. Existing native-direct/2 records and continuation routes are retained.

Native mixed field/electronic/LD/SBE persistence also uses CheckpointTarget/execute_checkpoint. Formats and native-direct/2 metadata are retained; callback failures stop subsequent slots, and completion-publication errors are agreed on all ranks.

The ordinary sample's `--direct-resume` defaults to protected restoration of field-local LD, SBE, native TDDFT/fast-PR HSE and Python clocks/currents after all groups prepare candidates. Old ABIs/direct-only targets are rejected before LOAD, without fallback. Explicit legacy behavior is `--unprotected-direct-restore`. See [CHECKPOINTING.en.md](CHECKPOINTING.en.md) for collective continuation/retry and failure exclusions.

`--memory-budget-gib 4` supplies an aggregate4GiB advisory budget. Live/candidate coexistence is reported and excess never stops the run. See [MEMORY.en.md](MEMORY.en.md).

## Mixed sample identity and restart

`python -m salmonpy.examples.mpi_mixed_native` preserves `--run-id`, logical point IDs, material IDs, material-definition and initial-state SHA256 digests, MPI ownership, and the LD cell-material mapping. The default `mixed-example` is a fixture label; provide an explicit run ID for each experiment. Output and temporary directories are excluded from definition fingerprints.

The same identity appears in `identity.json`, normal-run `result.json`, and the `point-response` archive. Identified direct checkpoints use `salmonpy.native-direct/3`; replay checkpoints use `salmonpy.native-replay/2`. Numerical payload formats and update order are unchanged. Direct restore checks identity and MPI ownership before loading numerical state; save checks agreement across all ranks.

Use `--legacy-checkpoint-identity` explicitly only for old unidentified checkpoints. This cannot bypass identity mismatches in identified checkpoints and does not disable protected restore. Custom hosts supply `ResponseSession` point records through `point_record` / `point_map`. These are caller declarations, not certification of physical accuracy.

## Native TDDFT rejection before update (limited scope)

Newly prepared and built standard embedded TDDFT exposes `sp_tddft_error_protocol() == 1`. Validation is collective on the attached material communicator:10 invalid/nonfinite input, stale A or missing input file;11 rank disagreement;12 unavailable lifecycle or repeated initialization;13 exhausted step capacity;14 pending restore candidate. These codes return unanimously before numerical mutation. Initial admission codes10/11 do not consume the lifetime; Python initialization can be retried with corrected input. Failed work directories remain for diagnosis.

Python maps only protocol-declared unanimous rejection to `InputError`. Old libraries, unknown codes or nonfinite results retain `RuntimeFailure`, with no inferred state guarantee. Every solver-group member must enter the same call. Invalid pointers, MPI/process failure, allocation errors, deep solver STOP and failures after numerical updates are outside this admission-only contract; see the limited deep-recovery contract below. Actual initialization still consumes the one-process lifetime. Fast HSE and native FDTD error contracts are unchanged.

This guarantee concerns an individual TDDFT state's pre-update rejection. A coupled session may already have advanced another component; it remains failed under the existing driver policy. This does not implement general coupled-step rollback.

## Expanded grid and parallel checks on a larger Mac

Keep the4GiB/4-core beginner recipe as its own verified profile. Ask for actual available resources and reserve headroom for the OS and other work. For the supplied64GB/18-core Mac, this campaign uses a48GiB advisory budget and at most15 MPI ranks with one thread each. These are campaign settings, not universal defaults or hard limits.

`mpi_mixed_native --field-shape 96 96 96` changes Maxwell grid counts: each axis must be even and >=40 and match its MPI decomposition. Default40³ and2-bohr spacing are unchanged. Enlarging the grid adds outer vacuum, so compare decompositions only at the same grid/input. Default material-point physical coordinates are preserved. Electronic GS grid/k sampling is unchanged: this is not a96³ TDDFT test. The3D fixture tests interfaces;1D remains recommended for production TDDFT coupling.

Field2+SBE1+TDDFT2 uses `mpiexec -n 5 ... --field-grid 1 2 1 --tddft-ranks 2`. Field12+SBE1+TDDFT2 uses `mpiexec -n 15 ... --field-grid 3 2 2 --field-shape 96 96 96 --tddft-ranks 2`. Replace `...` with the existing Python module invocation and BUILD,RUNTIME,GS,new-output arguments. Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.

Source-distribution tool `tools/test_mixed_scaling.py BUILD RUNTIME GS NEW_OUTPUT --rank-budget 18 --memory-budget-gib 48` uses the installed wheel and runs4/5/7/11/15 ranks sequentially twice on identical96³ inputs. A/E/j histories are compared with the first case. Different grid shapes are rejected. Supply MPI-vendor-specific flags with `--launcher-arg`; avoid other heavy work during timing.

`coupled_loop_seconds_max` excludes initialization, output and serial reference replay. Communication metrics include waits, not pure link bandwidth; native FDTD halo exchanges are not separately decomposed. Memory is the sum of per-rank peak RSS values, not simultaneous aggregate RSS or whole-system usage. Resource/timer measurement failures are distinct from calculation failures. These short repeats do not establish universal scalability or physical accuracy.

## Direct rollback after deep TDDFT numerical errors

New standard native TDDFT exposes `sp_tddft_recovery_protocol() == 1` and retains a per-step numerical payload in memory. Collective checks cover propagated wavefunctions(20), density(21), Hartree field(22), XC potential/energy(23), current/total energy(24), and a protected propagator error(25). Every material-group rank agrees before returning. Errors restore accepted orbitals, densities, fields, propagation history, XC state, counters and original overlap flags without propagation replay. Failed attempts add no ordinary RT output records.

A state returned by Python `NativeTDDFTBackend` raises `NumericalError` from `advance(A_old,A_new)` only for protocol-declared unanimous recovery. Its `A` and `steps` remain accepted; correct the cause and retry from that time. No automatic pulse weakening or dt changes are applied. Old libraries and unknown errors retain conservative `RuntimeFailure` handling without inferred recovery.

Supported: periodic, unpolarized, fixed-ion middlepoint PZ, automatic/fixed TB-mBJ, Libxc-PZ/PBE. SOC,+U,jellium, ionic motion/force output, Gram–Schmidt, projection/unfold and microscopic Maxwell coupling are outside the protected profile and return15 before updating. Ordinary batch calculations retain their existing route; external checkpoint encoding is unchanged. Per-step copies/checks add memory and runtime overhead: leave resource headroom. Budgets remain warning-only.

This is a TDDFT-group guarantee, not automatic rollback of an entire mixed FDTD step whose other components may have advanced. MPI/process death, allocation failures, uncontrolled third-party STOP, deep initialization STOP and output-device failures are excluded. One-rank fault injection exists only in a separately compiled `SALMONPY_RECOVERY_TESTING` build; production builds export no injection entry.

## Deep numerical recovery for fast-PR HSE

New libraries for `NativeFastHSEBackend` declare `sp_fast_hse_recovery_protocol() == 1`. The guarded profile is periodic, unpolarized, fixed-ion fast-PR HSE06 with Taylor4 / learned PT-CN and internal k-only MPI1/2. Conventional HSE is not a fallback.

Each step backs up electronic arrays and ACE, factor/PT-CN histories, exchange caches and counters. FFT plans, MPI communicators and other owning resources are not copied. Collective checks cover finite orbitals, density, Hartree, XC and endpoints, plus handled ACE operator errors, PT-CN convergence/Gram/retraction checks and history updates. Handled failures directly restore accepted numerical state without replaying propagation or exchange history. A, clocks, current and energy remain unchanged; failed attempts append no ordinary RT records. Diagnostic logs may remain.

Only a declared protocol and unanimous known error permit Python to raise `NumericalError` and retry with corrected input at the same accepted time. Old libraries and unknown errors retain conservative `RuntimeFailure` behavior. Neither dt nor drive amplitude changes automatically. Per-step backups/checks add memory and runtime cost; budgets remain advisory.

The guarantee covers the HSE material group, not automatic rollback of an entire mixed step containing FDTD and other materials. Deep initialization, allocation failure, process/MPI termination, uncontrolled external STOP and output-device failure remain excluded. Do not infer recovery after an unhandled failure.

Run `tools/test_fast_hse_deep_recovery_mpi.py BUILD GS NEW_OUTPUT --steps 32 --fault-step 13` with two ranks, OMP/BLAS1 and an external timeout. `--split` uses two independent two-rank groups. `--overflow-A 1e100` stresses real PT-CN numerical divergence; it is not a scientific input. `--faults` requires a separate explicit `SALMONPY_RECOVERY_TESTING` build; production exports no injection entry point. See `docs/fast-hse-deep-recovery-validation-2026-11-17.md`. These small synthetic cases verify implementation, not converged physical accuracy.

## Native FDTD numerical recovery

New standard host-MPI libraries declare `sp_eh_recovery_protocol() == 1`. `NativeEHField.advance` selects a checked entry and backs up field/PML, every field-local LD, deposited current and clock in memory. Collective gates check E after update (20), LD after update (21), H after update (22), and sampled fields/integrated A/ordinary observation arithmetic (23). Handled failures restore accepted state directly without propagation replay. A declared known recovery raises `NumericalError`, preserving Python A/E/clock/failed so the caller can correct input and retry. Failed attempts append no ordinary observation records; diagnostic/progress stdout may remain.

The guarded profile is the existing native E/H plane-wave path with local LD. TTM, dedicated linear response, ASE/ART and frequency-resolved observation are rejected before updating. Batch FDTD retains its old path. Old binaries without the checked entry and unknown failures remain conservative. Initialization, allocation, MPI/process termination, external STOP and output-device failure remain excluded. This protects the field group alone; it does not automatically roll back or retry a whole `CoupledStepper` / `ResponseSession` step. Budget for backup/check overhead; memory limits remain advisory.

Run `tools/test_eh_deep_recovery_mpi.py BUILD NEW_OUTPUT` with one/two ranks, OMP/BLAS1 and an external timeout. `--split` uses two independent two-rank groups. `--faults` requires a separate explicit `SALMONPY_RECOVERY_TESTING` build; production exports no injection entry. Overflow stress values are not scientific inputs.

## Explicit mixed-step protection

`CoupledStepper.step_protected(targets)` is an opt-in entry for complete backend-owned backup/rollback callbacks. See [STEPPING.en.md](STEPPING.en.md) and the six-rank `examples/mpi_protected_step.py`. Aggregate native transaction/output protection is not inferred from individual solver recovery.

Audited native aggregate protection uses the separate protocol and `--protect-steps` described in [NATIVE_STEPPING.en.md](NATIVE_STEPPING.en.md). Use `--hse-bundle` for distinct HSE initial data and `examples/mixed_tddft_hse_layers.json` for simultaneous TDDFT/HSE.

## Fast-PR HSE symmetry opt-in

NativeFastHSEBackend(...,allow_symmetry=True) is bounded to internal-k1/2,
Taylor/PT-CN with matching mode capability. The ordinary mixed sample explicitly selects --hse-symmetry
--hse-reference-bundle FULL_K_BUNDLE. Use --protect-steps/step_protected to roll
back the entire loop when InputError rejects a symmetry-breaking A. Direct/protected
resume and local two-LD/SBE/PBE-TDDFT/HSE mixtures are tested. Material/field/XC
compatibility and physical convergence remain separate; see the
[symmetry guide](FAST_HSE_SYMMETRY.en.md), including stable low-gradient evaluation and the original-case regression.

Reduced-k symmetry now supports two internal k ranks and same-partition direct resume with the separate native k-MPI capability. Checkpoint repartitioning remains unsupported; see the symmetry guide.

### Warnings versus execution errors

Memory/rank budgets, recommended scale and physical-accuracy advice are advisory: they do not block executable work. Inputs and timesteps are not automatically changed. Unreadable/corrupt files, incompatible dimensions or actual MPI partitions, and known unsupported conditions that stop native execution remain errors. Lack of validation alone is not a reason to refuse a calculation.

The fast-PR admission checks warn on finite occupations outside the expected range or a sum inconsistent with the electron count, without clipping or renormalizing them. Startup checks the orbital byte length without adding a full finiteness scan of potentially huge wavefunction files at every initialization. This does not guarantee detection of NaN orbitals before native startup. Legacy GS headers cannot authenticate physical GS identity.

SCF non-convergence is advisory. Finite saved states remain usable downstream; `GroundStateResult.converged=False` and saved residual/tolerance/iterations identify unconverged results.

Before fast-PR native startup, admission checks generated k meshes against retained symmetry, representative counts/weights, species pseudopotential references and GS header/occupation/orbital byte lengths. KY parsing accepts continuation records, unused trailing columns and implicit signed exponents, following native list-directed input. Other recognized PP formats receive reference/readability checks; native owns content parsing. Supported legacy GS records use native-endian4-byte markers/integers/logicals, binary64 occupations and one global real64/complex128 orbital stream. Arbitrary `file_kw`, unfolding and even Gamma endpoint meshes remain current fast-PR profile limitations. The original GS producer rank count need not match the new solver group; this differs from direct RT checkpoint partition compatibility.

Fast-PR HSE receives Cartesian three-component A/c, including x/y/z and circular polarization, without a fixed-z requirement. Reduced-k states check the actual vector against retained symmetry before each update; use a matching full-k initial state when no suitable group exists. The corrected core also supports general-cell symmetry initialization. See the [symmetry guide](FAST_HSE_SYMMETRY.en.md).

### Internal RT ownership (development)

Protected-step numerical buffers, output metadata and pending direct-restore
candidates now belong to their `tddft_context`, rather than a process singleton.
For internal non-symmetric fast-PR contexts sharing immutable model inputs and
one communicator, one context's rollback, staged commit/cancel or close does not
consume another's candidate. Different selected dt/capacity and three-component
waveforms can be interleaved. Reopened output streams acquire available units;
closing a rolled-back context cannot close a unit reused by another context.

The legacy Python create_state route retains one native lifetime/process.
The separate current C handle and Python open_model routes build on this internal
ownership; see SHARED_HSE.en.md. Different shared inputs, arbitrary communicators
and threaded entry remain separate work.
Shared symmetry/initialization/I/O/timer state and the synchronous serializer
allocation bound remain to be extracted. Do not copy whole contexts or owning
FFT/MPI/Libxc handles. Text-output rollback still restores saved file bytes on
the host before reopening; numeric checkpoints alone do not truncate logs.

Developers can explicitly run `tools/test_fast_hse_context_owners_mpi.py BUILD
FULL_K_FIXTURE NEW_OUTPUT`. It compiles an internal probe and compares isolated
and interleaved Taylor, learned Taylor and learned PT-CN executions, including
simultaneous protection, staged restore, cancellation, close and output-unit
reuse. Defaults are two independent solver groups with internal MPI1/2 and
OpenMP1/2. The fixtures test state isolation, not material convergence.

Reconnection failures retain the transaction metadata, so retry or close cannot act on another context’s reused unit.

Generation-checked same-model C contexts have a separate [native handle entry](NATIVE_HANDLES.en.md). The legacy Python path keeps its one-lifetime contract; [shared Python models](SHARED_HSE.en.md) use the separate handle route.


[Shared HSE multipoint placement and protected direct resume](SHARED_HSE.en.md#normal-mixed-mpi-sessions).

## Standard TDDFT initialization MPI retirement

A newly explicitly built standard runtime exposes
`sp_tddft_context_mpi_protocol()==1`. Communicators created during embedded TDDFT
initialization and its initial-current evaluation belong to that state. On close,
staged numeric restore candidates and protected buffers are discarded, request
caches and XC resources are finalized, and then its child communicators are
released. The host parent and unscoped attachment-owned children stay live until
their normal owner closes them. Checkpoint encoding and propagation order remain
unchanged. The **one native lifetime/process** restriction remains: retiring MPI
resources does not make native TDDFT reentrant or permit a second initialization.

For explicit ABI diagnostics, `int sp_tddft_resources(int values[4])` reports
tracked live children, registry slots, pending-restore-payload presence and live
TDDFT-state presence. Call on the main thread in the subgroup's normal order.
The small registry may retain vacant slots until detach; zero live children is
not a claim that all MPI/process memory is zero. These are experimental native
diagnostics, not a required method for custom Python backends.

PZ, Libxc-PZ/PBE and TB-mBJ automatic/fixed-c were compared against the preceding
runtime with two orbital or two real-space ranks. Currents match exactly; close
with a pending candidate returns tracked children/payload/state counts to zero
and leaves host collectives usable. Ordinary mixed full-k PBE/HSE/FDTD direct
restart and protected failed-restore retry remain covered. MPI/process failures,
initialization/allocation failures and different-input multi-context support are
separate work. Older runtimes retain their existing behavior; rebuild explicitly
when this retirement feature is required, rather than patching a fixed library.

MPI compiler wrappers must select an OpenMP-capable compiler. If a wrapper uses
a C compiler without OpenMP, select a supported compiler through that MPI's wrapper
configuration before building; the wrapper name alone does not identify its
underlying compiler. Record this choice with the build command.

The one-lifetime statements above describe the legacy Python/singleton route.
[Standard native handles](STANDARD_HANDLES.en.md) add a separate same-model C ABI;
they are not yet connected to a Python shared factory or ordinary multipoint plans.

For opt-in standard TDDFT shared models and multipoint groups, see [SHARED_TDDFT.en.md](SHARED_TDDFT.en.md).

For explicit sequential heterogeneous standard-model switching in one process, see [TDDFT_REGISTRY.en.md](TDDFT_REGISTRY.en.md).

New directory-scope3 runtimes avoid changing the host working directory during native E/H initialization, updates, checkpoint/restore and close. Relative files resolve against the field work directory; absolute paths remain unchanged, with scope cleanup on exceptions. Scope1/2 and older runtimes retain the previous directory-switching path. Numerical order and checkpoint formats are unchanged. Main-thread sequential calls and one native field lifetime per process remain required.

Standard native TDDFT shutdown releases its owned RT output streams and child communicators without closing a host-owned file on Fortran unit30. This is fixed in newly built runtimes; prepare a matching runtime explicitly rather than altering fixed binaries. This does not establish general I/O-failure recovery or threaded concurrency.

New runtimes keep the initialization-assigned MPI child owner active during standard TDDFT and fast-PR HSE propagation. Lazily created child communicators then retire with that state, preserving siblings and host MPI. Existing MPI decomposition, caller API and numerical checkpoint format are unchanged. This does not establish ownership of every resource, MPI-failure recovery or threaded concurrency.

New standard TDDFT singleton runtimes also isolate init/step/close timers from the host. With legacy timer protocol1, `state.timing_seconds` returns rank-local completed intervals, outside the numerical checkpoint. Older runtimes retain normal calculations; only this query requires the additive capability. Fast-HSE singleton timing has its corresponding legacy timer capability.

FDTD finite-value checks scan every array element through scalar IEEE inquiries. Collective rejection of NaN and positive/negative infinity, and state preservation on rejection, remain unchanged. No threshold or stopping policy was added. This revision is verified on GNU/macOS; Fujitsu compilation was not available.

New native FDTD close releases pending restore candidates and protected-step copies before closing the field. Early close still skips completion-only analysis such as final spectra. Numerical updates, checkpoint format and existing close status codes are unchanged. Optional read-only `sp_eh_resources(int flags[3])` reports field/protected-copy/pending-candidate allocation flags, not total memory or MPI resource retirement.

R2 native FDTD assigns child communicators created during initialization, updates and finalization to the field lifetime, releasing them after communication-cache retirement. Parent MPI and other owners remain intact. Optional read-only `sp_eh_mpi_resources(int values[3])` returns live connection children, registry slots and field-owned children. Slots may remain for reuse; live count determines release. Decomposition, E/H/LD update order and numerical checkpoint format are unchanged, as are the one-field-lifetime/process and MPI-failure-recovery limits.

Native FDTD initialization distinguishes admission from numerical setup. New runtimes collectively reject invalid point counts, empty/overlong/missing input paths, and successfully parsed unsupported theory/restart requests without consuming the native numerical lifetime. Parsing uses a local scratch workspace. An empty context is released if allocation or MPI-owner admission fails. Corrected native calls can retry on the same attachment before `eh_init`; the Python constructor attachment policy is unchanged. Malformed namelists and failures inside `eh_init` are not yet recoverable through this mechanism.

For new native FDTD runtimes, namelist read errors reported by the parser now return initialization failure collectively rather than finalizing MPI. Parser scratch files are closed and corrected input can retry before numerical setup. The ordinary standalone parser retains its previous failure behavior. This covers detected namelist read errors; earlier unit selection, coordinate parsing, semantic checks and numerical initialization still have separate failure paths.

Native FDTD now returns status1001 for an unsupported `unit_system` instead of
stopping MPI in the status-returning parser path. Python reports the accepted
spellings (`au`, `a.u.`, `A_eV_fs`). Parser scratch units close before rejection,
and no field state/child communicator is retained. Standalone parser callers
without an optional status retain their old stop behavior; valid unit conversion
and numerical update order are unchanged. This is not a general coordinate or
semantic-validation recovery contract, nor automatic job input repair.

The status-returning native FDTD parser also distinguishes coordinate-source
conflicts (1002/1003), reduced coordinates for an isolated system (1004), failed
coordinate-file open (1005) and failed row reads (1006). Root-rank read failures
are broadcast before coordinate-array collectives; external units and parser
scratch are closed. Python names the corresponding reason. These are existing
invalid-input conditions, not new physical restrictions. Valid coordinate
conversion/order and standalone no-status failure behavior are retained. This
scope does not cover allocation failure, every semantic check or arbitrary
native initialization failure. It is independent of the optional static checker.

Native E/H validates the y/n flags it actually uses (periodicity, restart,
shape/copy and field-observation/output flags) and returns1013 for invalid flags.
Uppercase Y/N normalization remains accepted. Electronic propagator, SCF and LCFO
choice lists are not applied to native Maxwell admission: values for those
inactive settings must not block an otherwise valid field calculation. This does
not mean Maxwell executes a named electronic method. Earlier common parsing and
active E/H numerical initialization still have their own checks; this is not a
promise to accept arbitrary unknown namelist keys or malformed physical inputs.
No-status standalone callers retain the existing global checks.

Custom Python/C current backends do not pass through the SALMON input parser and
are not restricted by its built-in method names. Their explicitly registered
backend owns parameter semantics; SALMONpy still checks the declared interface,
units/time contract and returned shapes/finite values at the connection. The
separate static checker treats unknown backend names/JSON parameters as
unverified, not unsupported, and never executes the supplied code.

Native E/H also skips electronic-only preprocessing during input reading: KY
pseudopotential/spin probing, rejection of conflicting EXX/legacy MLWF aliases,
EXX k-algorithm resolution and ACE warmup validation/arithmetic. Those inactive
settings and files must not prevent a field-only initialization. Input broadcasts
and field units/flags remain in place. This is scoped to status-returning
Maxwell reads; standalone no-status and electronic reads retain the existing
processing. Unknown namelist keys and malformed values still require appropriate
parser support; an ignored electronic method is not executed by the field solver.
