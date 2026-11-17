# MPI TDDFT coupled to Multiscale1D or native E/H FDTD (experimental)

[日本語](COUPLED.md) · Topic: `coupled` · Paired revision: `2026-11-14.1`

Split one Python MPI job into a field group and material groups, with **one material
and one macro point per material group**. Every member of a material group participates
in native TDDFT; k-point, orbital and real-space decomposition stays inside that group.
There is no subprocess or nested MPI launch. Native FDTD retains many local LD states.

The entries are `NativeTDDFTMaterialProvider`, `NativeTDDFTBackend`, `NativeEHField`
and `CoupledStepper`; the latter two live in `salmonpy.coupled`. Existing serial-worker
TDDFT and batch FDTD remain available with their own restrictions.

## Choosing dimension and response model

Three-dimensional fields and macro-point layouts are possible, but TDDFT keeps an
electronic state at every macro point. **Start with one-dimensional propagation or
layered structures for practical TDDFT work.** This refers to the electromagnetic/
macro layout, not reducing the microscopic DFT/TDDFT electronic problem to 1D.

**For three-dimensional calculations with electronic feedback, recommend SBE or a
user-written current solver.** Three-dimensional TDDFT is allowed, but choose it only
after assessing the required physics, independent state count, memory and runtime.
SBE/custom models must still satisfy their physical validity range, independent
point-state ownership and A/c-to-charge-current-density contract. This recommendation
does not imply an implemented SBE coupling adapter or automatic 1D reduction.
Standalone three-dimensional LD FDTD is a separate workflow.

The current 40³-grid/two-TDDFT-point sample is a connection regression, not a general
3D TDDFT recommendation or a transversely uniform 1D multilayer implementation.
Simply reducing transverse grid counts to one does not implement a 1D solver.

## Preparation and resources

First establish intent, total memory allowance and CPU concurrency. The standard
sample uses two field ranks plus two two-rank material groups: six ranks and two
electronic states. Start with one OpenMP/BLAS thread per rank. Memory estimation and
hard enforcement are not implemented; measure a small GS before scaling.

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
  synchronization. Post-native errors invalidate the state; no retry or rollback.
- Collectively close response, field, then Python context. Fortran does not finalize
  MPI. Native hangs/failures require launcher timeouts; the sample aborts the job on error.
- Restart, nonzero initial fields, impulse sources, rotations, unequal macro volumes,
  subcycling and TTM are not supported by this entry.

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

Every electronic grid cell gets its own state and MPI group; LD cells stay on field ranks. These are finite rectangular regions, not automatically homogenized infinite films. Increasing transverse area or thickness can greatly increase electronic ranks. The planner rejects excess ranks, unknown IDs, incompatible backends, absorbing-boundary overlap and electronic points in an LD constitutive neighborhood. It inserts no implicit gap. The rank budget is not a memory limit (`memory_bytes` is unknown). Initial-state preparation, MPI startup and physical convergence remain explicit.

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

Scope is periodic middlepoint propagation with the same environment/decomposition: PZ, TB-mBJ (automatic/fixed cval), Libxc-PZ/PBE. Libxc requires an explicitly built compatible native library. Generic custom-backend persistence, changed MPI partitions and other XC profiles require further implementation/validation; fast-PR HSE has a separate entry below. Abort the job on native restore failure; this is not atomic rollback of an existing live session. Allow storage and memory for full states.

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

This assigns field1,SBE1,HSE2; both LD materials stay local to the field rank. Use three separate jobs with the same arguments: continuous, `--direct-save-at 127`, then `--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT`. Rank-local HSE shards include group size/rank and k ranges; all shards are checked before committing numerical state. Fresh initialization recreates static exchange/FFT resources, not elapsed propagation. Restore requires the same binary, bundle, communicator topology and shared archive directory. Early invalid/truncated shards are rejected before numerical commit. A late I/O failure invalidates the contexts; abort the job rather than retrying a live session. No cross-group atomic rollback is promised.

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

Python maps only protocol-declared unanimous rejection to `InputError`. Old libraries, unknown codes or nonfinite results retain `RuntimeFailure`, with no inferred state guarantee. Every solver-group member must enter the same call. Invalid pointers, MPI/process failure, allocation errors, deep solver STOP and failures after numerical updates are excluded. Actual initialization still consumes the one-process lifetime. Fast HSE and native FDTD error contracts are unchanged.

This guarantee concerns an individual TDDFT state's pre-update rejection. A coupled session may already have advanced another component; it remains failed under the existing driver policy. This does not implement general coupled-step rollback.

## Expanded grid and parallel checks on a larger Mac

Keep the4GiB/4-core beginner recipe as its own verified profile. Ask for actual available resources and reserve headroom for the OS and other work. For the supplied64GB/18-core Mac, this campaign uses a48GiB advisory budget and at most15 MPI ranks with one thread each. These are campaign settings, not universal defaults or hard limits.

`mpi_mixed_native --field-shape 96 96 96` changes Maxwell grid counts: each axis must be even and >=40 and match its MPI decomposition. Default40³ and2-bohr spacing are unchanged. Enlarging the grid adds outer vacuum, so compare decompositions only at the same grid/input. Default material-point physical coordinates are preserved. Electronic GS grid/k sampling is unchanged: this is not a96³ TDDFT test. The3D fixture tests interfaces;1D remains recommended for production TDDFT coupling.

Field2+SBE1+TDDFT2 uses `mpiexec -n 5 ... --field-grid 1 2 1 --tddft-ranks 2`. Field12+SBE1+TDDFT2 uses `mpiexec -n 15 ... --field-grid 3 2 2 --field-shape 96 96 96 --tddft-ranks 2`. Replace `...` with the existing Python module invocation and BUILD,RUNTIME,GS,new-output arguments. Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.

Source-distribution tool `tools/test_mixed_scaling.py BUILD RUNTIME GS NEW_OUTPUT --rank-budget 18 --memory-budget-gib 48` uses the installed wheel and runs4/5/7/11/15 ranks sequentially twice on identical96³ inputs. A/E/j histories are compared with the first case. Different grid shapes are rejected. Supply MPI-vendor-specific flags with `--launcher-arg`; avoid other heavy work during timing.

`coupled_loop_seconds_max` excludes initialization, output and serial reference replay. Communication metrics include waits, not pure link bandwidth; native FDTD halo exchanges are not separately decomposed. Memory is the sum of per-rank peak RSS values, not simultaneous aggregate RSS or whole-system usage. Resource/timer measurement failures are distinct from calculation failures. These short repeats do not establish universal scalability or physical accuracy.
