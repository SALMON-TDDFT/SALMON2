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
  impulse sources, rotations, unequal macro volumes, subcycling remain unsupported.

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

The guarded profile is the existing native E/H plane-wave path with local LD. Dedicated linear response, ASE/ART and frequency-resolved observation are rejected before updating. Batch FDTD retains its old path. Old binaries without the checked entry and unknown failures remain conservative. Initialization, allocation, MPI/process termination, external STOP and output-device failure remain excluded. This protects the field group alone; it does not automatically roll back or retry a whole `CoupledStepper` / `ResponseSession` step. Budget for backup/check overhead; memory limits remain advisory.

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

Native FDTD initialization distinguishes admission from numerical setup. New runtimes collectively reject invalid point counts, empty/overlong/missing input paths, and successfully parsed unsupported theory/restart requests without consuming the native numerical lifetime. Parsing uses a local scratch workspace. An empty context is released if allocation or MPI-owner admission fails. Corrected native calls can retry on the same attachment before `eh_init`; the Python constructor attachment policy is unchanged. The detected parser errors and early initialization parameter errors described below are also recoverable. Other failures inside `eh_init` remain outside this mechanism.

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

Native FDTD now returns four existing early `eh_init` errors without terminating MPI: inconsistent grid specification (1101), periodicity/boundary conflict (1102), inconsistent time specification (1103), and CFL violation (1104). No new rejection criterion is added. The failed field context is discarded; if grid communication was prepared, its caches are released before its owned child communicators. A corrected native call can retry on the same attachment. Python constructor attachment ownership is unchanged.

MPI registry slots may remain for reuse: check live/field-owned counts, not zero slot capacity. Standalone callers omitting optional status keep their original diagnostics/termination. Recovery is limited to these early returns, not failures inside `eh_mpi_grid_sr`, later material/source/I/O stages, general allocation failure, OOM or MPI faults. No checker or calculation is automatically run on a management node.

The field parser route is selected explicitly with `read_input(..., status=..., native_host=.true.)`. Passing status alone does not skip ordinary preprocessing/checks. `native_host` defaults to false: returned errors and backend scope are separate choices.


## Spatially owned native FDTD current routing

With a matching runtime, CoupledStepper and PacketCoupledStepper query
`sp_eh_current_mask` once after field initialization. A point current is sent by
Scatterv to every field rank owning an edge touched by its transpose sample-map
deposit, rather than only the rank owning the point center. Boundary-crossing
points reach both sides; empty owners receive zero current entries. The three
components of an owned point travel together; this is not component-level compression.

The field leader retains globally ordered J for output and checkpoints. Other
field ranks retain the global array shape, but only owned rows are meaningful.
Use leader output for full-point analysis. A/E samples, edge weights, the
E→LD→H update order, collective failure decisions and restore guarantees are unchanged.

The new `sp_eh_step_owned` ABI accepts packed local point currents. The existing
`NativeEHField.advance(J)` still requires identical full-point J on every field
rank. Older runtimes without this ABI, Multiscale1D and standalone ResponseSession
keep replicated results. Build explicitly; importing or starting a calculation
does not rebuild a runtime automatically.

Checkpoint formats and MPI compatibility rules are unchanged. Full-row J saved
with the same compatible runtime and decomposition can resume by selecting the
owned rows. This does not bypass binary or MPI-layout compatibility checks.
The field manifest records `current_routing` and `current_rows_by_field_rank`;
explicit `communication_report()` rank records expose `received_current_rows`.
These are logical point counts, not measured network bytes.


## Coalesced outer step admission

Ordinary CoupledStepper combines completed field status and ResponseSession input
admission. That coalescing stage reduced successful-step outer agreements from 4 to 3.
The current standard entry further reduces them with status-bearing transport below.
At that stage PacketCoupledStepper decreased from 7 to 4 without H, or from 7 to 5 with H.
One agreement here is a hierarchical group→leaders→group decision, not a single
MPI call. Counts exclude setup, native validation/physical sums/halo exchange,
per-point subgroup checks and protected-step capture/rollback votes. They do not
state a wall-clock speedup.

The pre-field clock/usability decision and post-material decision remain. With H,
all groups agree on A/E sampling success before entering collective H sampling.
Field, sampling and packet-construction failures are rejected everywhere before
input distribution or material callbacks. Upstream errors travel as explicit
call arguments, not persistent session/process flags.

Standalone ResponseSession.step/step_fields signatures are unchanged. Custom
public response-step overrides and response objects without the merge hook keep
the previous separate agreement and dispatch order. Dispatch is negotiated once
across all groups at initialization: any rank with a custom entry selects the
legacy path everywhere. The coupling fixes response-entry dispatch when created.
Custom material backend contracts are unchanged.

Ordinary step is not newly transactional. Use the existing protected entry with
explicit StepTargets for whole-state preservation/retry, checking
original_preserved and cleanup_errors. Handled-exception tests do not establish
recovery from process loss, broken MPI or unsupported native failures.


## Status-bearing A/J transport

The standard coupled entry carries zero for success, or a reporting-rank code
for failure, in A/packet and J headers. Numeric A uses leader Alltoallv to
exchange source statuses and destination-owned A rows together; packets use
object alltoall. J and source status travel together through Gatherv. Field
subgroups receive their spatially owned J and result status in one Scatterv.

Only after every group knows failure from the headers do all participants enter
one explanatory-message Bcast, rooted at the selected reporting rank. No such
message Bcast occurs on normal steps. A failing rank never enters a Bcast while
healthy ranks proceed elsewhere unaware of it. Multiple reports select the
lowest reporting rank; that identifies a reporter, not necessarily the unique
physical cause of the failure.

ParallelContext duplicates a failure-only communicator once at creation and
frees it at close. It does not directly depend on MPI_COMM_WORLD or the caller's
input-communicator lifetime. Normal steps do not communicate on this domain;
import does not initialize MPI or create communicators.

Relative to the preceding implementation, separate hierarchical outer decisions
decrease A3→1, packets without H4→2, with H5→3. Actual Python leader collectives
also decrease5→4,6→5,7→6 respectively. Counts exclude native MPI and protected
transaction operations. The short post-J outcome notification remains necessary
for remote material groups to know whether to commit or restore; normal
communication is not zero.

Alltoallv status headers cost O(groups) per leader and O(groups²) in total. Fewer
calls do not establish multi-node speedup. Timing of status processing moves
into A_distribution/j_collection rather than a separate material decision.

Standalone ResponseSession and custom public response entries retain their prior
path. Unanimous setup selection, the pre-H vote, recipient-packet validation,
physical sums and protected capture/rollback decisions remain. Protected retry
guarantees and ordinary step's nontransactional semantics are unchanged. This
is not recovery from broken MPI or process loss.

### Retry after an undersized local MPI grid

New native E/H returns status1106 consistently on every field rank when a local
MPI subdomain violates the common grid initializer's existing stencil-width
minimum. This converts an existing STOP; it adds no new admissibility criterion,
partition rule or update change. Every rank already owns the complete domain
bounds, so rejection requires no extra agreement collective. It returns before
communication-cache creation and releases the failed state/owned MPI children.

For the local-grid case, unanimous1106 plus `sp_eh_init_recovery_protocol()==1` lets NativeEHField
close the rejected attachment and release its one-lifetime initialization
reservation. The caller keeps its MPI communicator and may construct again with
corrected settings and a new work directory. No such guarantee is inferred for
old libraries or unknown statuses. Resources/inputs are never adjusted
automatically: choose a layout/grid compatible with the job's allocation.
Standalone callers omitting status keep their existing diagnostic/termination.

Local GNU MPI1/2 tests include a domain too small on only the last rank. Two
low-level and two Python-constructor rejections permit successful initialization,
32 steps and direct checkpoint/restore. Live children return to0, caller MPI
survives, and normal A/E/H plus point output match the prior runtime exactly.
This does not cover all initialization/allocation failures, external STOP, OOM
or MPI/process failure.

### Recoverable explicit MPI layout/mapping errors

The same cold-initialization cleanup now covers status1105 (the existing Maxwell
k/orbital-factor requirement or a partition product inconsistent with the attached
solver group) and1107 (the active communicator mapper reaches its unimplemented
process_allocation branch). Recognized unanimous1105/1106/1107 and protocol1 let
the Python adapter detach and release the rejected lifetime reservation. Correct
the configuration and use a new work directory on the same caller communicator.
The existing startup status agreement is reused; no new collective is added.

Mapping names are not pre-whitelisted by Python or an outer validator. Returns
are placed in the mapper's actual error branches, preserving implemented modes,
Fujitsu mapping selection and caller-added native branches. Status-absent callers
retain existing diagnostics/termination. Automatic layout factorization failures,
Fujitsu topology errors, allocation and MPI/process failures remain separate
uncovered paths; this is not full initialization recovery or new support for an
unimplemented mapper. No scheduler/resource/input edits are automated.

GNU MPI1/2 regressions combine the existing admissions rather than duplicate the
suite:27 rejection/retry checks per layout, then normal32-step/direct-resume
verification with zero live children after cleanup. A/E/H, point output and
checkpoint bytes match the prior runtime. The existing orbital_sequential mode
is also explicitly exercised on MPI2; it remains usable.

### Recoverable material-ID CUBE input errors

Native E/H initialization returns1115 for unavailable, malformed, truncated,
nonfinite or grid-count-mismatched SALMON material-ID CUBE input. Reported tiny
reader staging-buffer allocation failure returns1116; this allocation path has
not been fault-injected and does not cover general OOM or solver allocations.
The reader owns its file unit and checks every value used in integer conversion,
including the final partial record. Existing streaming order and material-ID
conversion are preserved; this is not a general electronic-density CUBE reader.

For external shape input only, one field-group initialization agreement precedes
halo exchange, so a read failure on one rank returns on all ranks. No per-step
collective or generated-shape collective is added. Known unanimous statuses and
recovery protocol1 permit cold cleanup and retry with a new work directory on
the same caller communicator; status-absent standalone behavior is retained.

Local GNU MPI1/2 tests cover retry, caller MPI survival, direct resume and exact
agreement with prior generated-shape fields, point output and checkpoint bytes.
A host-owned unit600 survives rejected shape reads. Ordinary point headers,
point appends, final observation information, plane snapshots/integrals and the
material-ID CUBE writer now open a local newunit and close only that owned unit.
Filenames, record formats, ordering and numerical/communication operations are
unchanged. The generated-shape writer no longer closes the caller's requested
unit for an unimplemented output format. Existing helper signatures are retained.

GNU MPI1/2 tests keep a host-owned unit600 open throughout initialization,
32 steps, direct checkpoint restore and finalization. Field arrays, point files,
checkpoint bytes and generated shapes match prior controls; enabled plane output
also matches. This is file-unit ownership isolation, not output-error recovery:
write failures may still use existing Fortran error handling. Linear-response/spectral/ASE/ART units are addressed below; TTM has been removed;
legacy restart/helper ownership is addressed below. No claim covers arbitrary device failure,
external STOP, process failure or all initialization allocations.

### Legacy EM restart/helper file-unit ownership

The common_maxwell text/binary restart readers and writers now open their own
local newunit, including per-rank data and metadata files and the remapping read
path. txtfile_copy_em independently owns its input and output units instead of
borrowing700/701. Existing helper signatures retain ifn for source compatibility;
it no longer selects an actual file unit in these helpers. Paths, formats,
metadata, array order, copy trimming and communication/reconstruction operations
are unchanged. No additional startup or per-step collective is introduced.

The existing I/O fixture has an explicit --fdtd-restart-io mode: synthetic EM
arrays only, no SCF/GS data. GNU MPI1/2 tests preserve host units600/601/700/701
while exercising text values, real/complex1–5D binary arrays, single/all modes,
same-layout and cyclically exchanged target-domain reads, inactive payloads and
text copying. Control/occupied outputs and prior-runtime outputs match byte for
byte. This tests the helpers and does not certify every complete legacy batch
restart configuration or arbitrary rank counts. Normal LD+FDTD/direct-resume
field, point, plane, shape and checkpoint data also match prior controls.

Linear-response/spectral/ASE/ART units are addressed below; TTM has been removed.
This change isolates file-unit ownership; missing/corrupt restart data, IO device
failure and general allocation failure retain their existing handling. No general
IO-error recovery, Fujitsu execution or multi-node certification is claimed.

### Spectral output file-unit ownership

E/H finalization owns local newunit handles for linear-response trajectories and
spectra, frequency-resolved E/H/J planes, and ASE/ART spectra. The two simultaneous
windowed/unwindowed ASE/ART files use two local handles, not base-unit arithmetic.
The periodic linear-response response file is explicitly closed after its last
record; its earlier implicit close at process termination is no longer needed.
No module-wide file selector is added. Fourier transforms, normalizations, field
updates, loop/file ordering, names, record formats and MPI operations are retained.

The existing native I/O tool has --fdtd-spectrum-io, using actual small LD+FDTD
lifecycles with no electronic GS. GNU MPI1/2 controls cover isolated and periodic
impulse response, isolated non-unit-background raw response, Fourier planes, and
separate ASE and ART cases with three frequencies and both window variants. Caller-owned
units600/601/602 remain connected to their original files. Named spectral outputs
are closed before return; expected numeric outputs are nonempty and finite.
Unoccupied, occupied and previous-runtime outputs are compared byte for byte.
These are output/host-contract regressions in atomic units, not convergence or
material-property certification. Native E/H's existing supported-mode limits are
not lifted: the batch host used here enables spectral modes, not the embedded
step interface. No extra runtime collective, automatic build or startup run is
introduced. TTM has been removed; general IO-error recovery remains separate.

TTM is excluded from SALMONpy: unpublished implementation, automatic input-file activation, build targets, examples and utilities are removed from the maintained solver. A leftover `ttm.inp_ttm` file does not activate any calculation. Historical provenance records are retained.

### Recoverable automatic MPI factorization failure

Native E/H returns1117 when the existing automatic MPI distributor cannot assign
its rank count. The current automatic domain routine factors2/3/5; this is a limit
of that algorithm, not a blanket rank-count restriction. A supported explicit
spatial partition remains usable (tested with7 ranks). No Python/outer factor
whitelist is added and the distribution algorithm is unchanged. Correct the
partition and create a new work directory on the same caller communicator.

The internal Fortran helper has an optional status (0success,1unassigned prime
factor,2no allowed non-orthogonal k/orbital split). Without status it retains the
legacy diagnostics/STOP behavior. Matching Fortran modules must be used when
building callers; this is source compatibility, not stable Fortran binary ABI.
The C initialization signature and recovery protocol1 are unchanged. Python
releases a failed lifetime only for a recognized unanimous1117 with protocol1.
The existing cold cleanup/status agreement is reused; no new collective or
per-step work is added. This does not handle allocation, Fujitsu topology, MPI
failure or every specialized initialization stage.

Local GNU tests cover unchanged automatic choices on MPI1/2 and helper failure/
continuation on1/2/7. Native MPI7 repeats raw/public automatic-failure rejection,
keeps parent MPI alive, releases owned children and retries explicit decomposition;
32-step fields, point bytes and direct checkpoint/resume match the prior runtime.
This is local contract/regression evidence, not multi-node scaling certification.

### Retrying rejected LD tensor/subpixel input

The existing native readers already return1110 for missing/invalid diagonal or
reciprocal-tensor material sidecars and1113 for a missing/invalid reference
subpixel profile. For these cold errors the native initializer destroys the
candidate and releases its child resources. Python now recognizes these two
unanimous statuses with sp_eh_init_recovery_protocol1, detaches, and releases its
failed lifetime reservation. Correct/stage the material file and retry in a new
work directory on the same caller communicator. Error messages identify the
material input instead of reporting a generic native rejection.

No Fortran routine, input validation, coefficient arithmetic, file format,
communication or native binary changed in this correction. The existing audited
runtime is reused by recorded SHA256. Unknown/mixed statuses or an absent/wrong
recovery protocol do not enable retry. Handled coefficient/Bloch-mode returns are addressed below; arbitrary allocation,
IO-device, MPI or process failure remain outside the bounded guarantee.

The existing E/H recovery tool accepts --material-input-mode diagonal|tensor|subpixel
and --material-input-admission. GNU MPI1/2 repeat missing, invalid-header, truncated
and nonfinite-file rejection twice each through the public Python constructor;
candidate/attachment resources are released and the parent communicator remains
usable. Valid material initialization then completes32 steps and direct restore.
Compare against the pre-change Python package with the same native binary; this
is an input/lifecycle regression, not a new material-convergence certificate.

### Retrying handled coefficient/special-mode initialization errors

Python also recognizes existing native cold returns1111 (nonfinite diagonal/
tensor LD update coefficients),1112 (nonzero Bloch phase with the unsupported
off-diagonal tensor mode), and1114 (invalid subpixel update coefficients).
Only unanimous codes with initialization recovery protocol1 release the failed
lifetime. The native initializer already deallocates its candidate, releases its
send/receive cache before its owned MPI children, and returns; the correction is
Python's classification and diagnostic text, not a native numerical change.

Correct the input/mode and use a new work directory on the same communicator.
This does not implement Bloch/off-diagonal tensor evolution, silently set Bloch
phase to zero, relax finite checks or change coefficient formulas. Files already
written during rejected initialization remain in the rejected work directory;
the rejection is not a promise of filesystem rollback. Unknown/mixed statuses,
missing protocol, trapped floating-point exceptions and allocation/device/MPI/
process failures remain outside this bounded retry guarantee.

--coefficient-admission in the existing E/H tool uses material-input-mode
(diagonal,tensor,subpixel) to exercise1111/1112/1114. Finite but oversized test
rates provoke the existing coefficient guards; a native periodic test input
provokes the existing unsupported Bloch/tensor return. GNU MPI1/2 repeat each
failure and then run valid32-step/direct-resume controls. Reader rejections are
combined in the same fixture. The old SHA256-verified native binary is reused;
no solver rebuild, extra numerical calculation in Python or new collective.

### Retrying handled early grid/boundary/time errors

Native E/H already returns1101 for inconsistent cell/spacing/grid-count input,
1102 for incompatible periodic/boundary faces,1103 for an invalid duration/
step-count/time-step specification, and1104 when dt is not below its existing
CFL bound. Python now releases the failed initialization reservation for these
recognized unanimous returns only when initialization recovery protocol1 exists.
The native cold candidate/owned-child cleanup remains the same. Correct the
input and retry with a new work directory on the same caller communicator.

No new input restriction, automatic repair, CFL threshold change, native binary
change or collective is introduced. This aligns Python's error classification
with existing native behavior. Unknown/mixed statuses, missing protocol and
other unhandled failures remain conservative. Rejected directories are retained;
this does not promise rollback of staged or diagnostic files. No automatic job
resubmission or solver execution is added to preparation/checking workflows.

The existing --grid-init-retry fixture now includes these four public-constructor
rejections twice each. Combined with --geometry-admission, GNU MPI1/2 complete35
rejection/retry checks per run, verify parent MPI/child cleanup, then32 steps and
direct checkpoint restore. Field/point/shape/checkpoint bytes match earlier valid
controls. This is a local lifecycle regression, not all initialization/IO-error
recovery or cross-platform/physical-convergence certification.

### Core 3D allocation failure (1118)

Native E/H initialization uses `allocate(...,stat=...)` for the material-ID and
material scratch arrays and the core split E/H, coefficient and saved-E/H 3D
arrays. An allocation failure skips the remaining allocations in that phase.
The owning FDTD group agrees before shape/halo communication or downstream
coefficient setup and returns1118. Partially allocated arrays, communication
caches and owned child communicators are released. Python permits constructor
retry only for unanimous known1118 returns with the initialization recovery
protocol; unknown failures retain conservative lifetime handling.

Two phase agreements are added during initialization, within the field group.
There are no extra time-step communications. Memory-budget estimates still warn;
only an actual allocation failure ends this call.

The separate test build's `--allocation-admission` injects a STAT failure on one
rank after partial allocation and checks cleanup and retry. Production builds
contain no injection hook. LD pole, PML, observation and other auxiliary arrays
remain outside this bounded support. This does not guarantee recovery from OS
termination, physical memory exhaustion or MPI failure.

### LD state allocation failure (1119)

In addition to the core 3D arrays above, the initial `ld_allocate` now returns
STAT failures for its24 owned arrays:18 pole-current/history/polarization/
coefficient4D arrays and6 FDTD-current/old-E3D arrays. The caller owns partially
allocated state; native initialization cleanup releases it. One owning-field-
group agreement returns1119 on every rank if any rank failed. Python permits
retry only after unanimous known1119 and the initialization recovery protocol.
There are no extra time-step communications; memory-budget excess still warns.

Use `--ld-allocation-admission` only with a separately built testing runtime.
It fails one rank partway through pole and3D history allocations, then checks
cleanup and retry. Production includes no injection hook. LD update arithmetic
and order are unchanged. Material-list/coefficient scratch allocations, tensor
auxiliary arrays, subpixel-specific reallocation, PML and observation arrays
remain outside this bounded support, as do OS memory exhaustion/termination and
MPI failure.

### Material/boundary coefficient allocation failure (1120)

Native initialization detects STAT failures for three material vectors
(permittivity, permeability, conductivity) and six E/H boundary-coefficient
vectors. A partial allocation failure skips remaining allocations in that phase.
The material-vector gate shares the existing LD-state startup agreement.
This adds one field-group boundary agreement when LD is used (two startup
agreements without LD), before coefficient communication/vector initialization. Any failed rank returns1120
through the existing owned cleanup. Python constructor retry requires unanimous
known1120 and the initialization recovery protocol. No step collective is added;
memory-budget excess remains advisory.

The large PML3D coefficients are already covered by1118. PML coefficient setup
uses existing arrays and small automatic workspaces for the existing eight
layers. Arbitrary automatic/auxiliary allocation or OS memory failure is not
covered. PML formulas, layer count and update order are unchanged. The explicit
`--aux-allocation-admission` test build fails one rank after partial material or
boundary allocation; production contains no fault-injection hook.

### Optional point-observation I/O failures

Fortran checks open/write/close status for headers and time-series appends to
`obsN_at_point_rt.data`. A failure emits `SALMONPY_POINT_OUTPUT_WARNING` and
continues field/LD/electronic evolution. Warnings are limited to the first
failure per state and output-owning rank. Only successfully opened owned units
are closed; caller units remain untouched. No MPI communication is added.

The next observation retries the append. Once the destination is repaired and
a file with the required header exists, later rows can be written again. A
failed header creation requires preparing that header file. Missing or partially
written rows are not reconstructed or recomputed. Inspect warned output for
incompleteness/malformed rows before analysis.

This support is limited to point headers and appends. Plane integrals,info files,
restart data and other outputs are separate paths. It does not guarantee progress
through a filesystem-wide failure/full disk affecting required outputs. The
warning state is not part of the physical checkpoint format or numerical state.

### Optional plane-integral I/O failures

Fortran checks open/write/close for `obsN_xy/yz/xz_integral_rt.data` headers and
appends. Failure emits `SALMONPY_PLANE_INTEGRAL_WARNING` and continues evolution.
Integration,unit conversion and existing MPI sums are unchanged. A failed write
never skips subsequent plane or peer-rank reductions. No MPI communication is
added solely for this I/O handling.

The field state owns the bounded warning flag. Later observations retry after
the destination/header is repaired;missing/partial rows are not reconstructed.
Verify completeness before interpreting warned integral/energy-balance data.
The flag is excluded from the physical checkpoint format/numerical state,and
caller-owned file units remain untouched.

This support covers only integral headers/time-series appends,not snapshots,
frequency-resolved maps,info/restart files or all filesystem failures. Existing
Fortran source calls remain valid;`eh_save_plane_integral` accepts an optional
trailing warning flag. Without it,errors still warn,but cross-call warning
suppression is not retained.

### Optional observation-info failure during finalization

Fortran checks open/write/close for the animation/observation helper file
`obs0_info.data`. Failure emits `SALMONPY_OBSERVATION_INFO_WARNING` and proceeds
with normal state/MPI cleanup. A failed optional metadata write does not turn
completed physical evolution into a solver failure. Existing time series and
checkpoints are not removed;normal formatting remains unchanged. No MPI
communication or physical recomputation is added.

The info file can be missing/incomplete:check its presence/content before using
it for visualization or analysis. There is no automatic post-finalization retry.
This policy is limited to optional observation metadata;it does not silently
accept a failed required restart save or cover all output/filesystem failures.

### Optional plane-snapshot I/O failures

Fortran checks open/write/close for `obsN_<component>_<plane>_<step>.data`.
Failure emits `SALMONPY_PLANE_SNAPSHOT_WARNING` and continues other planes,
components and physical evolution. Existing plane aggregation,unit conversions,
MPI ordering,workspace release and normal formatting are preserved. Failed
writes never skip later reductions;no communication is added. The field owns
the bounded warning flag.

Later observations try their own time-stamped files normally. Missing/partial
files are not repaired automatically. After destination repair,the user can
advance from a checkpoint to the same time and write that snapshot again;the
restore itself adds no replay. Check snapshot availability/completeness before
visualization or analysis after a warning.

This scope covers time-resolved snapshot I/O,not staging-array allocation,
frequency-resolved maps,required restarts or filesystem-wide failures. Existing
Fortran source calls remain valid;`eh_save_plane` accepts an optional trailing
warning flag.

### Batch frequency-resolved plane-map I/O failures

Fortran checks open/write/close for
`obsN_eneM_<plane>_with_wf/without_wf.data`. Failure emits
`SALMONPY_FREQUENCY_MAP_WARNING` and continues remaining windows,planes,energies
and finalization. Fourier accumulation,unit conversion and MPI aggregation order
are unchanged;an output error never skips a later aggregation. No communication
or physical recomputation is added. Only successfully opened owned units are
closed;warnings are bounded within each finalization.

This support applies to the SALMON batch FDTD path. The current native E/H
coupling API's unsupported-diagnostic return15 is unchanged. Native protected
stepping and save/restore with these accumulators are not newly supported.
Missing/partial maps are not reconstructed or retried after finalization:check
energy/plane/window coverage before analysis. Accumulator/staging allocation,
ASE/ART and other output families,required restarts and filesystem-wide failures
remain outside this bounded support.

### Optional batch LR/ASE/ART spectra

Fortran checks open/write/close for LR trajectories/response tables and ASE/ART
spectral tables. Failure emits `SALMONPY_OPTIONAL_SPECTRUM_WARNING` and continues
remaining outputs and cleanup. ASE/ART keep independent units,status and names
for each window:one failed table never skips another window's integration,MPI
aggregation or output. Failed LR trajectory I/O does not skip Fourier processing
or the subsequent response table. Numerical formulas,units and normal formatting
are unchanged;no MPI communication or physical recomputation is added. Warnings
are bounded within each finalization.

This is the SALMON batch path;native E/H's unsupported-diagnostic return15 is
unchanged. Missing/partial spectra are not reconstructed or retried after
finalization. Check trajectory/energy/window coverage before analysis. This
optional-output policy does not silently accept failed required restart saves
or guarantee recovery from allocation/filesystem-wide failures.

### Required native checkpoint-save failures

Unlike optional diagnostics,`sp_eh_checkpoint(...,0)` returns nonzero on save
failure. It uses `status='new'` without replacing existing files and agrees
failure within the owning group. Detached payloads leave live field/LD state
and clocks unchanged. After catching the error,the host can repair the destination
or select a new one and retry from the same state. Never label a failed required
save as successful.

The low-level API does not provide atomic multi-file publication. New files
written by healthy ranks still belong to a failed set and are not an approved
checkpoint. Preserve committed data and quarantine failed-attempt files. Prefer
a new empty destination for retry;never delete/replace prior committed data to
reuse a path. Ranks must supply distinct names,and a higher layer publishes only
after every save and integrity check succeeds.

Same-configuration/topology normal states are verified with rank-local existing
file/directory collisions and missing parents,state retention,repaired saves and
direct restore. This validates existing recovery;no numerical routine or format
change was needed. Mid-stream full-disk,close-only,process/MPI failures and every
legacy batch `data_for_restart` writer remain outside the measured scope.

### Legacy batch EM restart-file saves

Legacy `data_for_restart` text/binary writers stage into an exclusively created adjacent `.salmonpy-pending` file rather than truncating the destination. Only successful open/write/close is followed by rename publication. Failed pending files remain as evidence. Before retry, move the affected pending file aside and retain existing good data.

The internal Fortran writers require `io_status`, initialized to zero by the caller at the start of each save phase. The first error is retained and skips subsequent local saves. Each initial-metadata/checkpoint phase makes one failure agreement within its owning MPI group. Required-save failure exits the batch with a nonzero status rather than reporting success. Ordinary non-checkpoint steps gain no agreement. This unpublished API has no omitted-status legacy fallback.

This is per-file staging. Binary payload publication waits until sidecar writing/close succeeds, but there is no atomic transaction across files/ranks or fsync durability guarantee. A failed phase can leave mixed generations; do not restart from that set. Retain a separately completed restart directory. This is distinct from native FDTD protected restoration using independent candidates.

### Batch restart-directory preparation

An existing `data_for_restart` directory returns success status zero; symlinks to valid directories are also accepted. Directory-creation failure replaces the previous post-creation barrier with one owning-group status agreement, followed by nonzero termination on all participating ranks rather than a root-only stop. A conflicting existing file or symlink target is left intact. This handles creation of a required destination; it is not a writability preflight or a guarantee that later writes succeed. Writers check actual I/O failures. Ordinary-step communication and numerical updates are unchanged.
