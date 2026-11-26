# Couple electronic calculations and light propagation from Python (MPI, experimental)

[日本語](COUPLED.md)

Feed the current produced by electrons back into FDTD to calculate the interaction
between light and matter. Python advances both calculations one step at a time.

Split one MPI job into processes responsible for the electromagnetic field and
processes responsible for electronic states. Each material group's k-point,
orbital and real-space parallelism stays within that group. The field group also
manages the LD states. No additional `mpiexec` is started during the calculation.

A macro point is a location on the electromagnetic grid where a material response is calculated. Start with one material and one macro point per material group. To use the same
material at several points, continue to model sharing; to handle different
materials in one group, continue to model registration. The conditions are listed below.

Use `salmonpy.coupled.NativeEHField` for the field, `NativeTDDFTMaterialProvider`
and `NativeTDDFTBackend` for the electrons, and `CoupledStepper` to connect them.
These are API names. Use `TDDFTResponseBackend` or ordinary FDTD execution when
each calculation should run in a separate process.

## Find the instructions for your task

| Goal | Read |
|---|---|
| Try an electronic calculation coupled to FDTD | [Preparation](#preparation-and-resources) → [Example](#sample) |
| Build layers from material IDs, order and thickness | [Layer setup](#assemble-ordered-layers-before-launching-mpi) |
| Use the same material at several points | [TDDFT](SHARED_TDDFT.en.md) or [fast HSE](SHARED_HSE.en.md) |
| Use different materials in one MPI group | [TDDFT model registration](TDDFT_REGISTRY.en.md) or [HSE model registration](SHARED_HSE.en.md) |
| Save/restart or retain the previous state after failure | [Checkpoints](CHECKPOINTING.en.md) or [failed time steps](NATIVE_STEPPING.en.md) |

## Choose how many calculation states to create

Choose the electronic-calculation entry by the number of states and material inputs.
The names below are Python methods.

| Goal | Entry | State handling |
|---|---|---|
| Start with one material and one point | `backend.create_state(...)` | This method allows one initialization per process |
| Apply different pulses or use several points for the same material | `backend.open_model(...)` → `model.new_state(...)` | Create independent states from shared material inputs |
| Use different material inputs in one MPI group | `backend.open_registry(...)` → `registry.open_model(...)` | Retain inputs per model and call them sequentially in the same group |

Operate the states sequentially from the main thread; this does not enable concurrent
calls from multiple threads. Check each dedicated guide and calculation library for
supported functionals, symmetry and MPI decompositions. Each field process can initialize
one FDTD field successfully. Retrying failed initialization is possible only for errors the library can recover from.

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
sample connects SBE to FDTD updates of E and H; dimension selection remains explicit.
Standalone three-dimensional LD FDTD is a separate workflow.

The current 40³-grid/two-TDDFT-point sample is a connection regression, not a general
3D TDDFT recommendation or a transversely uniform 1D multilayer implementation.
Simply reducing transverse grid counts to one does not implement a 1D solver.

## Preparation and resources

First establish intent, total memory allowance and CPU concurrency. The standard
sample uses two field ranks plus two two-rank material groups: six ranks and two
electronic states. Start with one OpenMP/BLAS thread per rank. The mixed sample offers uncertain memory estimates and measurements; excess only warns. Measure a small GS before scaling.

Explicitly build a **shared library against the MPI managed by the calling program**. Portable serial runtimes
and the communication-only bridge cannot substitute for it. Nothing builds on import.
The current build route targets macOS/Linux and GNU-compatible compilers.

```sh
python tools/prepare_salmon_runtime.py NEW_SOURCE
python tools/build_embedded_tddft.py NEW_SOURCE NEW_BUILD --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=mpifort \
  --cmake-arg=-DCMAKE_C_COMPILER=mpicc \
  --cmake-arg=-DBLA_VENDOR=OpenBLAS \
  --cmake-arg=-DOPENMP_FLAGS=-fopenmp \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --cmake-arg='-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -ffree-line-length-none -fallow-argument-mismatch'
```

Resolve compiler/BLAS locations for the calling program. The C wrapper must accept the selected
OpenMP flags too. The tool enables MPI/PIC, audits sources and objects for forbidden
eigensolvers, and records library SHA256, source hashes and commands in `build.json`.
Tools also ship under `salmonpy/_runtime_source/tools`; locate them with
`importlib.resources.files('salmonpy').joinpath('_runtime_source')`.

mpi4py and the library must use the same MPI. This calling program-specific library is not a
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

Providers and response calculation implementations in `MaterialBinding` remain replaceable user objects.
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

The factory checks material-site A(0) before creating calculation states. A nonzero
value is rejected with its point/value: increase the left vacuum cells or delay
the incident pulse. Incident A in vacuum may be nonzero. Do not silently zero
material A(0) or substitute an initially illuminated GS. This first connection
supports only zero initial material A and one solver initialization per process.

The update remains previous j → Multiscale1D field step → new A → material step.
Initial feedback is zero. Both clocks are at n·dt (offset 0); the E/H half-step
adapter does not apply. `step()` returns point-ordered `A_over_c` and `J_charge`
on the field rank. Keep the original field state to read its `grid_A` and
`grid_E`; the latter retains its centered-difference time (n−1)·dt.
Each returned frame declares `E_time`: Multiscale1D's centered E is at
(n−1)·dt while A is at n·dt; FDTD updates of E and H reports E at n·dt. Do not align E and A
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
measurements explicitly. Single-calling program timings do not establish multi-node scaling.

Acceptance compares full-grid A and each material point's j against two-point
SALMON Weyl multiscale output, alongside serial-worker and rank-count parity.
A synthetic GS checks the interface/update order, not physical convergence.

## Field times and grid positions

The new scheme is `salmonpy.eh-centered-current/1`. It is a different discretization
from the Weyl multiscale path; bitwise equivalence between those schemes is not
claimed. FDTD updates of E and H ordering and incident-source evaluation times are retained.

| Quantity | Field clock |
|---|---|
| E after update n; solver saved E/H (H is time-averaged) | n·dt |
| New A and J; raw updated H | (n+1/2)·dt |
| Previous J consumed by E update | (n−1/2)·dt |

Electronic local time zero corresponds to field time dt/2. Initial E, H, A and
feedback J are zero. A nonzero computed initial GS current does not replace that
initial feedback. The sequence is previous J → FDTD updates of E and H update → A update →
TDDFT J update. `A_new = A_old − dt * S(E_n)`; crystal drive is always A/c.

This table is declared in `salmonpy.time/1` form (`NativeEHField.timing`: j in at
offset ½ with lag 1 and t=0 unused, A out at offset ½). `CoupledStepper` gathers
every point's matter declaration over MPI, places it on the electronic clock
(field clock + ½, `ELECTRONIC_CLOCK_OFFSET`) and rejects mismatches before the
first step; `stepper.timing` and `response.json` record both. See
[FIELDS.en.md](FIELDS.en.md).


With integer k counting dt/2, old A/psi/J/raw H have k=2n-1, E(n) and centered H(n) have k=2n, and new A/psi/J/raw H have k=2n+1. FDTD consumes previous J at2n-1; midpoint TDDFT uses the mean of old/new A. Since A_new-A_old=-dt*S(E_n), both use the same midpoint2n.

Crystal TDDFT uses the same midpoint A for the kinetic operator and nonlocal projector phase during propagation, then switches both to endpoint A before current/energy evaluation. See [A and observation step numbers](SHARED_TDDFT.en.md#crystal-tddft-a-and-observation-step-numbers). The prescribed-A propagation is unchanged; observations previously evaluated at inconsistent A values are corrected.

Current molecular one-way E/H input uses the same electronic clock correspondence. A two-frame transfer does not mean two half-step propagations. See the [molecular step table](ISOLATED_TDDFT.en.md#half-steps-and-update-numbers).

S matches solver observation collocation: Ex averages edges `(i,j,k)` and
`(i−1,j,k)`; Ey/Ez use the corresponding previous y/z edge. Deposit J using S transpose,
half to each corresponding edge. Each macro point represents one uniform grid-cell
volume, so no extra volume multiplier applies. Each field rank adds only to edges it owns.

Coupled points stay at least ten cells from boundaries and clear of the incident
source plane/support and nonvacuum constitutive edges. The surrounding 3³ background
cells must be vacuum. LD may occupy other regions.


## Advance light more often than electrons

`CoupledStepper` is the equal-step interface. Use
`FieldSubcycledStepper(response, field, field_steps=m)` for m field updates per
electronic update. Set electronic step delta_t=response.dt and field step
Delta_t=field.dt with the explicit integer ratio delta_t=m*Delta_t. Field-local
LD advances at each field step; external electronic states advance at delta_t.

Each cycle advances the field m times with the **last accepted current held**,
then sends endpoint A, advances the electrons once and receives new j. There is
no future-current interpolation, prediction or implicit iteration. The hold error
is part of the accuracy study. A/j are exchanged once per electronic update;
no extra application-wide communication is added per light update. All ranks
call step once per electronic update. Field-only before_field_step/after_field_step
callbacks can observe each actual light update.

For light index l, E is at l*Delta_t and raw H/A at(l+1/2)*Delta_t. Electronic
local zero corresponds to field Delta_t/2. New electronic state/current is at
n*delta_t+Delta_t/2; the electronic midpoint is(n-1/2)*delta_t+Delta_t/2.
**With unequal steps that midpoint need not coincide with a stored E time.**
The kinetic operator and projector phases still use the same midpoint A, and
observations use the same endpoint A. Current midpoint A is the mean of the two
electronic endpoints; the complete fine-field history is not passed to TDDFT.
Retain field_step/electronic_step and timing's two steps/time origin.

This interface supports periodic A-driven matter. Molecular E/H driving,
noninteger ratios, coarser light steps, whole-cycle protected restore and disk
restart are not validated here. Do not continue a failed composite cycle or infer
the equal-step interface's restoration guarantees for this new schedule.

## Choose TDDFT + FDTD timesteps

Set the electronic timestep δt and field timestep Δt separately. `FieldSubcycledStepper` requires δt=mΔt: advance the field m times, then update the electrons once. The field step must satisfy the CFL condition for the Maxwell grid; the ratio alone does not determine a suitable step.

Use the same GS, grids, k points, pulse and final physical time for comparisons. First fix the electronic step and halve the field step; then vary the electronic step separately. Compare currents at matching physical times. Evaluate absorption with **Σ Δt V J·E**, using the actual current supplied to FDTD and E averaged before/after each field update. Internal TDDFT energy gain is a separate auxiliary diagnostic.

A5% relative difference in current and absorption is an initial screening target. Tighten it for the intended application and check electron-number loss too. Field-step convergence does not establish convergence in electronic step, spatial grids or k sampling.

### Si comparison example

This example uses electronic dt0.04 au, fourth-order Taylor propagation, a16³ electronic grid,2×2×2 k points and a40³ Maxwell grid with20 bohr spacing. One FDTD point carries Si/PZ response. All cases reach20 atomic time units.

| Field dt (au) | Field updates / electronic update | Absorption (Ha, including macro-point volume) |
|---|---:|---:|
| 0.01 | 4 | 7.27652472e-5 |
| 0.005 | 8 | 7.13166003e-5 |
| 0.0025 | 16 | 7.05918506e-5 |

Adjacent relative current RMS differences are0.0231% and0.0114%. Relative absorption differences are2.03% and1.03%, normalized by the finer result. This short, weakly driven example meets the5% screening target; it does not prescribe steps for other materials, strong excitation or long runs.

## Stored fields, outputs and limits

- One calculation state lifetime per process. No reinitialization, multiple states in
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
  MPI. solver hangs/failures require launcher timeouts; the sample aborts the job on error.
- Same-topology protected restart is documented below. Nonzero initial fields,
  impulse sources, rotations and unequal macro volumes remain unsupported.
  Use the separate explicit FieldSubcycledStepper for multiple light updates.

`response.json` stores point-ordered A/c, sampled E, previous/new charge currents
and both clocks. `salmonpy-coupled-field.json` records mapping, MPI decomposition,
library and completion status. Material working copies and solver output files remain available.
Select large spatial outputs and sampling intervals in the RT input. These JSON
files are not complete restart checkpoints.

Local regressions check original vacuum/LD FDTD parity, current injection at MPI
boundaries in all three directions, one/two field-rank parity and TDDFT A-history
replay. Physical grid/time/GS convergence and production scaling remain separate work.

## Different materials and time steps

Material states in one ResponseSession share an electronic exchange interval,
regardless of their material/grid/XC. Default CoupledStepper also requires that
field step and rejects mismatched dt. Explicitly choose FieldSubcycledStepper
only for finer light steps with electronic_dt=m*field_dt. Finer electrons and
noninteger ratios are not supported by this interface. No intermediate-A,
time-averaged or future-current prediction is applied automatically. The physical
electronic origin is field_dt/2, not electronic_dt/2.

## Toward mixed LD, SBE and TDDFT

`mpi_mixed_multiscale` couples LD and a synthetic SBE point to one 1D field:
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_mixed_multiscale NEW_OUTPUT
```
One field rank and two material ranks check nonzero responses and exact agreement with an independent serial schedule. Output includes point IDs10/20, material keys, A/c, charge current j and communication timing. Endpoint LD derives E from the A difference; it is not the FDTD updates of E and H semi-implicit LD scheme.

With four ranks, `--runtime RUNTIME --bundle GS_BUNDLE --input rt.inp` prepares a third managed serial TDDFT worker point. Four-rank,100-step three-way execution using the Si starter GS (PZ,16³ grid,2×2×2 k mesh) matches an independent serial calculation exactly. Supply a compatible GS/runtime explicitly. This validates the managed serial worker, not internally MPI-parallel TDDFT or converged optical properties. Evidence: `docs/mixed-three-material-validation-2026-11-16.md`.

The final topology keeps many LD states local to FDTD ranks while electronic groups own SBE/TDDFT states. Do not copy the demonstration's one-LD-cell-per-rank layout into large production jobs. This is not real-material convergence or validated 3D mixed coupling.

## Field-local solver LD with SBE and internally parallel TDDFT

`mpi_mixed_native` uses one field rank holding two LD regions, one SBE rank and two ranks for one TDDFT point. Four outer ranks exercise orbital decomposition inside the TDDFT group.
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT
```
Use three outer ranks with `--tddft-ranks 1`. `--without-ld` replaces LD regions with vacuum. BUILD_JSON is an explicit `tools/build_embedded_tddft.py` build for the MPI managed by the calling program; SERIAL_RUNTIME provides an independent response replay. Nothing builds or downloads implicitly.

LD retains FDTD updates of E and H coefficients, auxiliary variables and update order. Electronic points retain A/c drive and charge current j. Point10 is SBE,20 TDDFT; constitutive IDs1/2 belong to the field material grid. Electronic neighborhoods do not overlap LD. Field half-step clocks are distinct from electronic local time.

This tiny40³/2-electronic-point 3D fixture is a connection regression; keep recommending 1D propagation for practical TDDFT. It does not establish large-scale LD performance, multiple distributed field ranks, whole-session restart or physical convergence. solver failures are not recoverable; the sample aborts its MPI job.

Validated:128 steps,1/2 internal TDDFT rank field/current comparison, serial response replay, vacuum-LD contrast, original LD-only parity and installed-wheel execution. Scope/errors: `docs/native-mixed-material-validation-2026-11-16.md`.

### Parallelize the field group too

Two field ranks along y, one SBE rank and one TDDFT rank:
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT --field-grid 1 2 1 --tddft-ranks 1 --compare REFERENCE_OUTPUT
```
REFERENCE_OUTPUT is the same sample with one field and one TDDFT rank (three outer ranks). Omit `--compare` to run independently. Comparison requires matching steps,dt,material configuration and point order, then checks A/c,E,j histories. Required outer ranks: `product(field-grid)+1+tddft-ranks`.

Two-way x/y/z field decompositions each matched the one-field-rank result exactly over128 steps. The y/z layouts split both LD materials across both field blocks. Simultaneous two-field/two-TDDFT-group ranks require five ranks and were not tested within this four-core campaign. This is not a large-scale performance or memory-cap guarantee. Evidence: `docs/distributed-field-mixed-validation-2026-11-16.md`.

## Assemble ordered layers before launching MPI

`salmonpy.layers.assemble_native_layers` accepts a `salmonpy.native-layers/1` layout and an explicit mapping from material IDs to LD models or `MaterialBinding` objects. Use existing configured material loaders/bindings for physical parameters and initial states; the layout references their IDs. Custom response calculation implementations use the same `MaterialBinding` contract.

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

Every electronic grid cell gets its own state and MPI group; LD cells stay on field ranks. These are finite rectangular regions, not automatically homogenized infinite films. Increasing transverse area or thickness can greatly increase electronic ranks. The planner warns on an exceeded rank budget and preserves every point and the required rank count. Execution still requires that many MPI ranks. It rejects unknown IDs, incompatible calculation implementations, absorbing-boundary overlap and electronic points in an LD constitutive neighborhood. It inserts no implicit gap. The rank budget is not a memory limit (`memory_bytes` is unknown). Initial-state preparation, MPI startup and physical convergence remain explicit.

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

To repeat SBE instead, change the condition above to `sbe-point`. Both layouts were checked over128 steps by independently replaying each point's distinct A/c history and comparing currents. This does not test multiple TDDFT called directly from Python states within one process.

### Replay restart for the coupled sample

`mpi_mixed_native --save-at 64` stops at step64 and writes input fingerprints and A/E/current histories to `OUTPUT/checkpoint/`. Use the same total steps, layout, library, GS and MPI topology with a new output directory and `--resume PREVIOUS_OUTPUT/checkpoint`. Add `--compare CONTINUOUS_OUTPUT` to compare complete histories.

This reconstructs field, LD and electronic histories by rerunning their original schedule from initial conditions. Recorded A/E/current must match exactly before continuation. It requires the original assets and repeats the elapsed computation; it is not a direct internal-state dump. Scope is the supplied sample in the same environment/topology. General coupled-state archives remain unimplemented; see the direct sample restart below.

### Direct internal-state restart

Use `mpi_mixed_native --direct-save-at 64` to stop at64 and save rank-local state in `OUTPUT/direct-checkpoint/`. Start a new output/job with the same inputs, library, GS and MPI topology and `--direct-resume PREVIOUS_OUTPUT/direct-checkpoint`. Keep the original total128 steps. No propagation to the saved time is repeated. `--compare CONTINUOUS_OUTPUT` compares histories; `result.json` reports `restored_step` and `native_steps_this_job`.

Saved state includes solver E/H/PML, LD P/J/Q and old fields, electronic orbital buffers/density/potentials/Hartree/RT histories/A/current, SBE density matrix, and Python integrated A/lagged current/clocks. Initialization rebuilds static coefficients/resources before loading dynamic arrays. A completed file inventory with SHA256 checks is required; changed inputs/library/topology are rejected.

Scope is periodic middlepoint propagation with the same environment/decomposition: PZ, TB-mBJ (automatic/fixed cval), Libxc-PZ/PBE. Libxc requires an explicitly built compatible Fortran/C library. Generic custom-calculation implementation persistence, changed MPI partitions and other XC configurations require further implementation/validation; HSE has a separate entry below. Protected candidate restore is now the default and preserves originals on handled staging failures; the explicit unprotected route has no rollback guarantee. Allow storage and memory for full states.

Direct checkpoint format v2 validates XC names, derivative settings and input cval while retaining freshly initialized Libxc handles. Preserve runtime/checkpoint hashes. The sample RUNTIME argument also accepts an explicitly built/audited serial SALMON executable for reference validation; do not use an MPI binary as this child process. Use matching GS and RT functional inputs for each XC.

#### Internally parallel TB-mBJ/Libxc-PBE direct restart

Use matching-XC GS/RT bundles and four outer ranks: field1,SBE1,TDDFT orbital ranks2. Supply `--tddft-ranks 2 --rank-budget 4 --layers LAYOUT_JSON`. Run a continuous reference, a save job with `--direct-save-at 63`, and a fresh output/job with `--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT`. Keep total steps128.

For TB-mBJ (automatic coefficient/fixed cval1.2) and Libxc-PBE, the65 resumed steps produced complete A/c,E,current histories bit-identical to continuous128. This validates direct same-topology two-level MPI restart. Changing decomposition remains unsupported. Field2+SBE1+TDDFT2 needs five ranks and remains outside this four-core campaign.

Checksum verification on direct restart is partitioned across ranks. Every rank must read the same completed manifest and all checks must succeed before state restoration. Payloads must reside in one shared directory visible with identical contents to all ranks; different node-local copies are outside this contract. `result.json` records per-rank file/byte counts in `checkpoint_validation`.

## HSE in outer MPI / FDTD coupling

`FastHSEMaterialProvider` / `FastHSEResponseBackend` connect a dedicated serial HSE worker to the existing A/c→charge-current j contract. FDTD keeps its solver MPI entry, LD stays on field ranks, and each HSE point's single-rank group owns an independent child process. This is outer MPI coupling, not internally parallel HSE. calling program cwd is unchanged; import starts no MPI, build or download.

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

Explicitly link the worker from an audited HSE build into a new external directory. Conventional HSE is not an engine option.

```sh
python tools/build_fast_hse_worker.py FAST_PR_BUILD_JSON NEW_WORKER_DIRECTORY
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_mixed_native FIELD_BUILD_JSON UNUSED_REFERENCE HSE_BUNDLE NEW_OUTPUT --fast-hse-build NEW_WORKER_DIRECTORY/build.json --layers PATH_TO_mixed_hse_layers.json --dt .001 --steps 256 --rank-budget 3
```

The packaged layout contains two LD materials,SBE,HSE. `UNUSED_REFERENCE` is an unused positional argument for this HSE-only example and may be a nonexistent placeholder; independent reference validation uses another identical worker. For two field ranks use `mpiexec -n 4 ... --field-grid 1 2 1 --rank-budget 4`. Each HSE group requires one rank. The child uses that group's CPU while its MPI caller waits. Budget extra CPUs explicitly before increasing OpenMP threads.

Supply an original GS/RT bundle with tddft_response,hse06 and an explicit fast propagator (hse_taylor4/hse_taylor4_full/hse_ptcn). Only selected dt/step capacity are overridden; other physics inputs remain. Use serial internal parallel settings and the fast engine's dense-k/fixed-occupation/full-support restrictions for learned modes. The zero-initial-field sample requires initial HSE current<=1e-12 au and subtracts no background current. The preparation tool generates a synthetic continuation fixture, not a converged material GS.

`--direct-save-at 127` captures field/PML/LD/SBE and worker internal state at one epoch. A new job uses `--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT`. LOAD restores accepted step/current, with saved A/coupling histories checked. Only129 remaining steps are propagated; no elapsed127-step reconstruction. Fresh static initialization and independent serial validation replay remain separate. Same worker,bundle,input,MPI topology and shared archive directory are required.

POSIX macOS/Linux worker uses /dev/stdin; measured here on Apple Silicon. Orbital/grid-parallel HSE, changed-topology restore,GPU,generic calculation implementation archives and continuation after solver failure remain unsupported.

Worker inputs must use atomic units (unit_system=au), no internal pulses: ae_shape1=none or explicit zero impulse, ae_shape2=none, trans_longi=tr. Double driving by an internal pulse and external FDTD A is rejected. State time is accepted step×dt.

## HSE with internal k-point MPI

`salmonpy.fast_hse_mpi.NativeFastHSEBackend` runs HSE collectively on an explicitly supplied material communicator. The calling program initializes MPI and retains ownership of finalization. Import does not load MPI or compile. Each group owns one point and one solver initialization per process, on the main thread. The call to the Fortran/C library temporarily changes the process working directory: do not use other threads that access relative paths.

Prepare a new explicit-routing source with `tools/prepare_hse_rt_state.py NEW_SOURCE --explicit-routing`, then explicitly build it with `tools/build_fast_hse_mpi.py PREPARED_SOURCE NEW_OUTPUT --jobs 2 --cmake-arg=...`. Use matching MPI managed by the calling program Fortran/C compilers, OpenBLAS/LAPACK, FFTW and optional Libxc paths as in the runtime build guide. This build requires MPI/HSE and rejects GPU, Poisson-FFTW, ScaLAPACK and EigenExa configurations. The validation record binds source and adapter hashes. Only the pinned HSE configuration is an engine option.

```python
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend
from salmonpy.material_response import MaterialBinding

binding = MaterialBinding("hse-point",
    FastHSEMaterialProvider(bundle_directory, "rt.inp", max_steps=256),
    NativeFastHSEBackend(mpi_build_json, existing_shared_work_root))
# Assign two ranks to this point with assemble_native_layers / ParallelPlan.
```

Supported decomposition is k-only, one or two ranks, `nproc_ob=1,nproc_rgrid=1,1,1`. The calculation implementation overrides `nproc_k` to the communicator size and overrides dt/step capacity explicitly. Use enough k points for the group. Other driving/fixed-occupation/full-support restrictions above remain. The sample requires a serial worker validation record too, for independent reference validation after the coupled run.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native FIELD_BUILD_JSON UNUSED_REFERENCE HSE_BUNDLE NEW_OUTPUT --fast-hse-build SERIAL_WORKER_BUILD_JSON --fast-hse-mpi-build MPI_BUILD_JSON --hse-ranks 2 --layers PATH_TO_mixed_hse_layers.json --dt .001 --steps 256 --rank-budget 4
```

This assigns field1,SBE1,HSE2; both LD materials stay local to the field rank. Use three separate jobs with the same arguments: continuous, `--direct-save-at 127`, then `--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT`. Rank-local HSE shards include group size/rank and k ranges; all shards are checked before committing numerical state. Fresh initialization recreates static exchange/FFT resources, not elapsed propagation. Restore requires the same binary, bundle, communicator topology and shared archive directory. Early invalid/truncated shards are rejected before numerical commit. Default protected restore preserves original states on handled staging failures. Explicit unprotected restore and unhandled solver failures have no retry guarantee.

Measured on Apple Silicon with a synthetic 8³ microscopic grid and 4×1×1 k mesh: Taylor/learned PT-CN, field1+SBE1+HSE2, 256 steps at dt=.001, direct split127. All A/c,E,j history arrays matched bit-for-bit after129 resumed calls. Against independent serial HSE response, maximum HSE current differences were4.07e-16 and1.71e-16 au. These are connection/restart regressions, not material convergence or scaling certification. Evidence: `docs/fast-hse-kmpi-validation-2026-11-16.md`.

During layer planning, calculation implementations may declare `supported_group_sizes=(1,2)`: a nonempty tuple of unique positive integers. Unsupported rank counts are rejected before MPI or initial-state preparation. Fast HSE declares `(1,2)`. user-written calculation implementations without this optional declaration retain existing execution/decomposition checks. This is not certification of numerical/restart capabilities and ResponseSession also validates it and agrees declarations across all groups before preparing materials.

The direct sample writes manifest schema `salmonpy.native-direct/2` with persistence format,save/restore modes,group rank/point ownership and field/electronic k-orbital-real-space decompositions. All ranks compare it with the current candidates before any LOAD, including layout changes with the same total rank count. Static resources initialize first; preflight performs no propagation,snapshot,restore or LOAD. Legacy manifest/1 lacks this metadata and is rejected by the new route; use the original compatible reader for old archives. solver payload formats are unchanged. Custom direct-state/1 callbacks require explicit checkpoint_format,checkpoint_modes (tuple of save/restore) and checkpoint_decomposition (positive-integer tuple whose product equals group size). Declarations are caller contracts, not certification of implementations or atomic rollback. A snapshot method alone does not establish direct restart support; the supplied helper accepts only built-in LD/SBE snapshot routes, not replay-based states.

Capability and MPI compatibility checks use the common `salmonpy.checkpointing` API. See [CHECKPOINTING.en.md](CHECKPOINTING.en.md) for user-written calculation implementations and serial hosts. Existing solver-direct/2 records and continuation routes are retained.

solver mixed field/electronic/LD/SBE persistence also uses CheckpointTarget/execute_checkpoint. Formats and solver-direct/2 metadata are retained; callback failures stop subsequent slots, and completion-publication errors are agreed on all ranks.

The ordinary sample's `--direct-resume` defaults to protected restoration of field-local LD, SBE, solver TDDFT/HSE and Python clocks/currents after all groups prepare candidates. Old ABIs/direct-only targets are rejected before LOAD, without fallback. The ordinary example has no unprotected-restore option. See [CHECKPOINTING.en.md](CHECKPOINTING.en.md) for collective continuation/retry and failure exclusions.

`--memory-budget-gib 4` supplies an aggregate4GiB advisory budget. Live/candidate coexistence is reported and excess never stops the run. See [MEMORY.en.md](MEMORY.en.md).

## Mixed sample identity and restart

`python -m salmonpy.examples.mpi_mixed_native` preserves `--run-id`, logical point IDs, material IDs, material-definition and initial-state SHA256 digests, MPI ownership, and the LD cell-material mapping. The default `mixed-example` is a fixture label; provide an explicit run ID for each experiment. Output and temporary directories are excluded from definition fingerprints.

The same identity appears in `identity.json`, normal-run `result.json`, and the `point-response` archive. Identified direct checkpoints use `salmonpy.native-direct/3`; replay checkpoints use `salmonpy.native-replay/2`. Numerical payload formats and update order are unchanged. Direct restore checks identity and MPI ownership before loading numerical state; save checks agreement across all ranks.

Old arguments, checkpoint schemas and missing-field defaults are not supported. Use the current API and checkpoint writer. The mixed example always records material and initial-state identities.

## TDDFT called directly from Python rejection before update (limited scope)

Newly prepared and built standard TDDFT called in the Python process exposes `sp_tddft_error_protocol() == 1`. Validation is collective on the attached material communicator:10 invalid/nonfinite input, stale A or missing input file;11 rank disagreement;12 unavailable lifecycle or repeated initialization;13 exhausted step capacity;14 pending separate state prepared for restoration. These codes return unanimously before numerical mutation. Initial input checking codes10/11 do not consume the lifetime; Python initialization can be retried with corrected input. Failed work directories remain for diagnosis.

Python maps only protocol-declared unanimous rejection to `InputError`. Old libraries, unknown codes or nonfinite results retain `RuntimeFailure`, with no inferred state guarantee. Every solver-group member must enter the same call. Invalid pointers, MPI/process failure, allocation errors, deep solver STOP and failures after numerical updates are outside this input checking-only contract; see the limited deep-recovery contract below. Actual initialization still consumes the one-process lifetime. Fast HSE and FDTD called directly from Python error contracts are unchanged.

This guarantee concerns an individual TDDFT state's pre-update rejection. A coupled session may already have advanced another component; it remains failed under the existing driver policy. This does not implement general coupled-step rollback.

## Expanded grid and parallel checks on a larger Mac

Keep the4GiB/4-core beginner recipe as its own verified configuration. Ask for actual available resources and reserve headroom for the OS and other work. For the supplied64GB/18-core Mac, this campaign uses a48GiB advisory budget and at most15 MPI ranks with one thread each. These are campaign settings, not universal defaults or hard limits.

`mpi_mixed_native --field-shape 96 96 96` changes Maxwell grid counts: each axis must be even and >=40 and match its MPI decomposition. Default40³ and2-bohr spacing are unchanged. Enlarging the grid adds outer vacuum, so compare decompositions only at the same grid/input. Default material-point physical coordinates are preserved. Electronic GS grid/k sampling is unchanged: this is not a96³ TDDFT test. The3D fixture tests interfaces;1D remains recommended for production TDDFT coupling.

Field2+SBE1+TDDFT2 uses `mpiexec -n 5 ... --field-grid 1 2 1 --tddft-ranks 2`. Field12+SBE1+TDDFT2 uses `mpiexec -n 15 ... --field-grid 3 2 2 --field-shape 96 96 96 --tddft-ranks 2`. Replace `...` with the existing Python module invocation and BUILD,RUNTIME,GS,new-output arguments. Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.

Source-distribution tool `tools/test_mixed_scaling.py BUILD RUNTIME GS NEW_OUTPUT --rank-budget 18 --memory-budget-gib 48` uses the installed wheel and runs4/5/7/11/15 ranks sequentially twice on identical96³ inputs. A/E/j histories are compared with the first case. Different grid shapes are rejected. Supply MPI-vendor-specific flags with `--launcher-arg`; avoid other heavy work during timing.

`coupled_loop_seconds_max` excludes initialization, output and serial reference replay. Communication metrics include waits, not pure link bandwidth; FDTD called directly from Python halo exchanges are not separately decomposed. Memory is the sum of per-rank peak RSS values, not simultaneous aggregate RSS or whole-system usage. Resource/timer measurement failures are distinct from calculation failures. These short repeats do not establish universal scalability or physical accuracy.

## Direct rollback after deep TDDFT numerical errors

New standard TDDFT called directly from Python exposes `sp_tddft_recovery_protocol() == 1` and retains a per-step numerical payload in memory. Collective checks cover propagated wavefunctions(20), density(21), Hartree field(22), XC potential/energy(23), current/total energy(24), and a protected propagator error(25). Every material-group rank agrees before returning. Errors restore accepted orbitals, densities, fields, propagation history, XC state, counters and original overlap flags without propagation replay. Failed attempts add no ordinary RT output records.

A state returned by Python `NativeTDDFTBackend` raises `NumericalError` from `advance(A_old,A_new)` only for protocol-declared unanimous recovery. Its `A` and `steps` remain accepted; correct the cause and retry from that time. No automatic pulse weakening or dt changes are applied. Old libraries and unknown errors retain conservative `RuntimeFailure` handling without inferred recovery.

Supported: periodic, unpolarized, fixed-ion middlepoint PZ, automatic/fixed TB-mBJ, Libxc-PZ/PBE. SOC,+U,jellium, ionic motion/force output, Gram–Schmidt, projection/unfold and microscopic Maxwell coupling are outside the protected configuration and return15 before updating. Ordinary batch calculations retain their existing route; external checkpoint encoding is unchanged. Per-step copies/checks add memory and runtime overhead: leave resource headroom. Budgets remain warning-only.

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

## Explicit mixed-step protection

`CoupledStepper.step_protected(targets)` is an opt-in entry for complete calculation implementation-owned backup/rollback callbacks. See [STEPPING.en.md](STEPPING.en.md) and the six-rank `examples/mpi_protected_step.py`. Aggregate solver transaction/output protection is not inferred from individual solver recovery.

Audited solver aggregate protection uses the separate protocol and `--protect-steps` described in [NATIVE_STEPPING.en.md](NATIVE_STEPPING.en.md). Use `--hse-bundle` for distinct HSE initial data and `examples/mixed_tddft_hse_layers.json` for simultaneous TDDFT/HSE.

## HSE symmetry opt-in

NativeFastHSEBackend(...,allow_symmetry=True) is bounded to internal-k1/2,
Taylor/PT-CN with matching mode capability. The ordinary mixed sample explicitly selects --hse-symmetry
--hse-reference-bundle FULL_K_BUNDLE. Use --protect-steps/step_protected to roll
back the entire loop when InputError rejects a symmetry-breaking A. Direct/protected
resume and local two-LD/SBE/PBE-TDDFT/HSE mixtures are tested. Material/field/XC
compatibility and physical convergence remain separate; see the
[symmetry guide](FAST_HSE_SYMMETRY.en.md), including stable low-gradient evaluation and the original-case regression.

Reduced-k symmetry now supports two internal k ranks and same-partition direct resume with the separate solver k-MPI capability. Checkpoint repartitioning remains unsupported; see the symmetry guide.

### Warnings versus execution errors

Memory/rank budgets, recommended scale and physical-accuracy advice are advisory: they do not block executable work. Inputs and timesteps are not automatically changed. Unreadable/corrupt files, incompatible dimensions or actual MPI partitions, and known unsupported conditions that stop solver execution remain errors. Lack of validation alone is not a reason to refuse a calculation.

The HSE input checking checks warn on finite occupations outside the expected range or a sum inconsistent with the electron count, without clipping or renormalizing them. Startup checks the orbital byte length without adding a full finiteness scan of potentially huge wavefunction files at every initialization. This does not guarantee detection of NaN orbitals before solver startup. Legacy GS headers cannot authenticate physical GS identity.

SCF non-convergence is advisory. Finite saved states remain usable downstream; `GroundStateResult.converged=False` and saved residual/tolerance/iterations identify unconverged results.

Before HSE solver startup, input checking checks generated k meshes against retained symmetry, representative counts/weights, species pseudopotential references and GS header/occupation/orbital byte lengths. KY parsing accepts continuation records, unused trailing columns and implicit signed exponents, following solver list-directed input. Other recognized PP formats receive reference/readability checks; solver owns content parsing. Supported legacy GS records use solver-endian4-byte markers/integers/logicals, binary64 occupations and one global real64/complex128 orbital stream. Arbitrary `file_kw`, unfolding and even Gamma endpoint meshes remain current HSE configuration limitations. The original GS producer rank count need not match the new solver group; this differs from direct RT checkpoint partition compatibility.

HSE receives Cartesian three-component A/c, including x/y/z and circular polarization, without a fixed-z requirement. Reduced-k states check the actual vector against retained symmetry before each update; use a matching full-k initial state when no suitable group exists. The corrected core also supports general-cell symmetry initialization. See the [symmetry guide](FAST_HSE_SYMMETRY.en.md).

### Internal RT ownership (development)

Protected-step numerical buffers, output metadata and pending direct-restore
candidates now belong to their `tddft_context`, rather than a process singleton.
For internal non-symmetric HSE contexts sharing immutable model inputs and
one communicator, one context's rollback, staged commit/cancel or close does not
consume another's candidate. Different selected dt/capacity and three-component
waveforms can be interleaved. Reopened output streams acquire available units;
closing a rolled-back context cannot close a unit reused by another context.

The legacy Python create_state route retains one solver initialization per process.
The separate current C handle and Python open_model routes build on this internal
ownership; see SHARED_HSE.en.md. For different model inputs in one MPI group, use [registered HSE models](SHARED_HSE.en.md). Switching to arbitrary communicators or concurrent calls from multiple threads is outside this interface.
Shared symmetry/initialization/I/O/timer state and the synchronous serializer
allocation bound remain to be extracted. Do not copy whole contexts or owning
FFT/MPI/Libxc handles. Text-output rollback still restores saved file bytes on
the calling program before reopening; numeric checkpoints alone do not truncate logs.

Developers can explicitly run `tools/test_fast_hse_context_owners_mpi.py BUILD
FULL_K_FIXTURE NEW_OUTPUT`. It compiles an internal probe and compares isolated
and interleaved Taylor, learned Taylor and learned PT-CN executions, including
simultaneous protection, staged restore, cancellation, close and output-unit
reuse. Defaults are two independent solver groups with internal MPI1/2 and
OpenMP1/2. The fixtures test state isolation, not material convergence.

Reconnection failures retain the transaction metadata, so retry or close cannot act on another context’s reused unit.

Generation-checked same-model C contexts have a separate [solver handle entry](NATIVE_HANDLES.en.md). The legacy Python path keeps its one-lifetime contract; [shared Python models](SHARED_HSE.en.md) use the separate handle route.


[Shared HSE multipoint placement and protected direct resume](SHARED_HSE.en.md#normal-mixed-mpi-sessions).

<details>
<summary>Internal resources, communication and validation records (for developers)</summary>

The material below supports building and extending calculation libraries. It also
records specific implementation stages and tested configurations. For ordinary
execution, use the instructions above; for shared or registered models, follow
their dedicated guides.

## Standard TDDFT initialization MPI retirement

A newly explicitly built standard runtime exposes
`sp_tddft_context_mpi_protocol()==1`. Communicators created during TDDFT called in the Python process
initialization and its initial-current evaluation belong to that state. On close,
staged numeric separate states prepared for restoration and protected buffers are discarded, request
caches and XC resources are finalized, and then its child communicators are
released. The calling program parent and unscoped attachment-owned children stay live until
their normal owner closes them. Checkpoint encoding and propagation order remain
unchanged. The **one solver initialization per process** restriction remains: retiring MPI
resources does not make TDDFT called directly from Python reentrant or permit a second initialization.

For explicit ABI diagnostics, `int sp_tddft_resources(int values[4])` reports
tracked live children, registry slots, pending-restore-payload presence and live
TDDFT-state presence. Call on the main thread in the subgroup's normal order.
The small registry may retain vacant slots until detach; zero live children is
not a claim that all MPI/process memory is zero. These are experimental solver
diagnostics, not a required method for custom Python calculation implementations.

PZ, Libxc-PZ/PBE and TB-mBJ automatic/fixed-c were compared against the preceding
runtime with two orbital or two real-space ranks. Currents match exactly; close
with a pending candidate returns tracked children/payload/state counts to zero
and leaves calling program collectives usable. Ordinary mixed full-k PBE/HSE/FDTD direct
restart and protected failed-restore retry remain covered. MPI/process failures,
initialization/allocation failures and different-input multi-context support are
separate work. Older runtimes retain their existing behavior; rebuild explicitly
when this retirement feature is required, rather than patching a fixed library.

MPI compiler wrappers must select an OpenMP-capable compiler. If a wrapper uses
a C compiler without OpenMP, select a supported compiler through that MPI's wrapper
configuration before building; the wrapper name alone does not identify its
underlying compiler. Record this choice with the build command.

The one-lifetime statements above describe the legacy Python/singleton route.
[Standard solver handles](STANDARD_HANDLES.en.md) add a separate same-model C ABI;
they are connected to the [Python shared-model and multipoint interfaces](SHARED_TDDFT.en.md).

For opt-in standard TDDFT shared models and multipoint groups, see [SHARED_TDDFT.en.md](SHARED_TDDFT.en.md).

For explicit sequential heterogeneous standard-model switching in one process, see [TDDFT_REGISTRY.en.md](TDDFT_REGISTRY.en.md).

New directory-scope3 runtimes avoid changing the calling program working directory during FDTD updates of E and H initialization, updates, checkpoint/restore and close. Relative files resolve against the field work directory; absolute paths remain unchanged, with scope cleanup on exceptions. Scope1/2 and older runtimes retain the previous directory-switching path. Numerical order and checkpoint formats are unchanged. Main-thread sequential calls and one solver field lifetime per process remain required.

Standard TDDFT called directly from Python shutdown releases its owned RT output streams and child communicators without closing a calling program-owned file on Fortran unit30. This is fixed in newly built runtimes; prepare a matching runtime explicitly rather than altering fixed binaries. This does not establish general I/O-failure recovery or threaded concurrency.

New runtimes keep the initialization-assigned MPI child owner active during standard TDDFT and HSE propagation. Lazily created child communicators then retire with that state, preserving siblings and MPI managed by the calling program. Existing MPI decomposition, caller API and numerical checkpoint format are unchanged. This does not establish ownership of every resource, MPI-failure recovery or threaded concurrency.

New standard TDDFT singleton runtimes also isolate init/step/close timers from the calling program. With legacy timer protocol1, `state.timing_seconds` returns rank-local completed intervals, outside the numerical checkpoint. Older runtimes retain normal calculations; only this query requires the additive capability. Fast-HSE singleton timing has its corresponding legacy timer capability.

FDTD finite-value checks scan every array element through scalar IEEE inquiries. Collective rejection of NaN and positive/negative infinity, and state preservation on rejection, remain unchanged. No threshold or stopping policy was added. This revision is verified on GNU/macOS; Fujitsu compilation was not available.

New FDTD called directly from Python close releases pending separate states prepared for restoration and protected-step copies before closing the field. Early close still skips completion-only analysis such as final spectra. Numerical updates, checkpoint format and existing close status codes are unchanged. Optional read-only `sp_eh_resources(int flags[3])` reports field/protected-copy/pending-candidate allocation flags, not total memory or MPI resource retirement.

R2 FDTD called directly from Python assigns child communicators created during initialization, updates and finalization to the field lifetime, releasing them after communication-cache retirement. Parent MPI and other owners remain intact. Optional read-only `sp_eh_mpi_resources(int values[3])` returns live connection children, registry slots and field-owned children. Slots may remain for reuse; live count determines release. Decomposition, E/H/LD update order and numerical checkpoint format are unchanged, as are the one-field-lifetime/process and MPI-failure-recovery limits.

FDTD called directly from Python initialization distinguishes input checking from numerical setup. New runtimes collectively reject invalid point counts, empty/overlong/missing input paths, and successfully parsed unsupported theory/restart requests without consuming the solver numerical lifetime. Parsing uses a local scratch workspace. An empty context is released if allocation or MPI-owner input checking fails. Corrected calls to the Fortran/C library can retry on the same attachment before `eh_init`; the Python constructor attachment policy is unchanged. The detected parser errors and early initialization parameter errors described below are also recoverable. Other failures inside `eh_init` remain outside this mechanism.

For new FDTD called directly from Python runtimes, namelist read errors reported by the parser now return initialization failure collectively rather than finalizing MPI. Parser scratch files are closed and corrected input can retry before numerical setup. The ordinary standalone parser retains its previous failure behavior. This covers detected namelist read errors; earlier unit selection, coordinate parsing, semantic checks and numerical initialization still have separate failure paths.

FDTD called directly from Python now returns status1001 for an unsupported `unit_system` instead of
stopping MPI in the status-returning parser path. Python reports the accepted
spellings (`au`, `a.u.`, `A_eV_fs`). Parser scratch units close before rejection,
and no field state/child communicator is retained. Standalone parser callers
without an optional status retain their old stop behavior; valid unit conversion
and numerical update order are unchanged. This is not a general coordinate or
semantic-validation recovery contract, nor automatic job input repair.

The status-returning FDTD called directly from Python parser also distinguishes coordinate-source
conflicts (1002/1003), reduced coordinates for an isolated system (1004), failed
coordinate-file open (1005) and failed row reads (1006). Root-rank read failures
are broadcast before coordinate-array collectives; external units and parser
scratch are closed. Python names the corresponding reason. These are existing
invalid-input conditions, not new physical restrictions. Valid coordinate
conversion/order and standalone no-status failure behavior are retained. This
scope does not cover allocation failure, every semantic check or arbitrary
solver initialization failure. It is independent of the optional static checker.

FDTD updates of E and H validates the y/n flags it actually uses (periodicity, restart,
shape/copy and field-observation/output flags) and returns1013 for invalid flags.
Uppercase Y/N normalization remains accepted. Electronic propagator, SCF and LCFO
choice lists are not applied to solver Maxwell input checking: values for those
inactive settings must not block an otherwise valid field calculation. This does
not mean Maxwell executes a named electronic method. Earlier common parsing and
active E/H numerical initialization still have their own checks; this is not a
promise to accept arbitrary unknown namelist keys or malformed physical inputs.
No-status standalone callers retain the existing global checks.

Custom Python/C current calculation implementations do not pass through the SALMON input parser and
are not restricted by its built-in method names. Their explicitly registered
calculation implementation owns parameter semantics; SALMONpy still checks the declared interface,
units/time contract and returned shapes/finite values at the connection. The
separate static checker treats unknown calculation implementation names/JSON parameters as
unverified, not unsupported, and never executes the supplied code.

FDTD updates of E and H also skips electronic-only preprocessing during input reading: KY
pseudopotential/spin probing, rejection of conflicting EXX/legacy MLWF aliases,
EXX k-algorithm resolution and ACE warmup validation/arithmetic. Those inactive
settings and files must not prevent a field-only initialization. Input broadcasts
and field units/flags remain in place. This is scoped to status-returning
Maxwell reads; standalone no-status and electronic reads retain the existing
processing. Unknown namelist keys and malformed values still require appropriate
parser support; an ignored electronic method is not executed by the field solver.

FDTD called directly from Python now returns four existing early `eh_init` errors without terminating MPI: inconsistent grid specification (1101), periodicity/boundary conflict (1102), inconsistent time specification (1103), and CFL violation (1104). No new rejection criterion is added. The failed field context is discarded; if grid communication was prepared, its caches are released before its owned child communicators. A corrected call to the Fortran/C library can retry on the same attachment. Python constructor attachment ownership is unchanged.

MPI registry slots may remain for reuse: check live/field-owned counts, not zero slot capacity. Standalone callers omitting optional status keep their original diagnostics/termination. Recovery is limited to these early returns, not failures inside `eh_mpi_grid_sr`, later material/source/I/O stages, general allocation failure, OOM or MPI faults. No checker or calculation is automatically run on a management node.

The field parser route is selected explicitly with `read_input(..., status=..., native_host=.true.)`. Passing status alone does not skip ordinary preprocessing/checks. `native_host` defaults to false: returned errors and calculation implementation scope are separate choices.


## Spatially owned FDTD called directly from Python current routing

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


## Coalesced outer step input checking

Ordinary CoupledStepper combines completed field status and ResponseSession input
input checking. That coalescing stage reduced successful-step outer agreements from 4 to 3.
The current standard entry further reduces them with status-bearing transport below.
At that stage PacketCoupledStepper decreased from 7 to 4 without H, or from 7 to 5 with H.
One agreement here is a hierarchical group→leaders→group decision, not a single
MPI call. Counts exclude setup, solver validation/physical sums/halo exchange,
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
Custom material calculation implementation contracts are unchanged.

Ordinary step is not newly transactional. Use the existing protected entry with
explicit StepTargets for whole-state preservation/retry, checking
original_preserved and cleanup_errors. Handled-exception tests do not establish
recovery from process loss, broken MPI or unsupported solver failures.


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
also decrease5→4,6→5,7→6 respectively. Counts exclude solver MPI and protected
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

New FDTD updates of E and H returns status1106 consistently on every field rank when a local
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
Fujitsu mapping selection and caller-added solver branches. Status-absent callers
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

FDTD updates of E and H initialization returns1115 for unavailable, malformed, truncated,
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
A calling program-owned unit600 survives rejected shape reads. Ordinary point headers,
point appends, final observation information, plane snapshots/integrals and the
material-ID CUBE writer now open a local newunit and close only that owned unit.
Filenames, record formats, ordering and numerical/communication operations are
unchanged. The generated-shape writer no longer closes the caller's requested
unit for an unimplemented output format. Existing helper signatures are retained.

GNU MPI1/2 tests keep a calling program-owned unit600 open throughout initialization,
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
arrays only, no SCF/GS data. GNU MPI1/2 tests preserve calling program units600/601/700/701
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

The existing solver I/O tool has --fdtd-spectrum-io, using actual small LD+FDTD
lifecycles with no electronic GS. GNU MPI1/2 controls cover isolated and periodic
impulse response, isolated non-unit-background raw response, Fourier planes, and
separate ASE and ART cases with three frequencies and both window variants. Caller-owned
units600/601/602 remain connected to their original files. Named spectral outputs
are closed before return; expected numeric outputs are nonempty and finite.
Unoccupied, occupied and previous-runtime outputs are compared byte for byte.
These are output/calling program-contract regressions in atomic units, not convergence or
material-property certification. FDTD updates of E and H's existing supported-mode limits are
not lifted: the batch calling program used here enables spectral modes, not the in-process
step interface. No extra runtime collective, automatic build or startup run is
introduced. TTM has been removed; general IO-error recovery remains separate.

TTM is excluded from SALMONpy: unpublished implementation, automatic input-file activation, build targets, examples and utilities are removed from the maintained solver. A leftover `ttm.inp_ttm` file does not activate any calculation. Historical provenance records are retained.

### Recoverable automatic MPI factorization failure

FDTD updates of E and H returns1117 when the existing automatic MPI distributor cannot assign
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
continuation on1/2/7. solver MPI7 repeats raw/public automatic-failure rejection,
keeps parent MPI alive, releases owned children and retries explicit decomposition;
32-step fields, point bytes and direct checkpoint/resume match the prior runtime.
This is local contract/regression evidence, not multi-node scaling certification.

### Retrying rejected LD tensor/subpixel input

The existing solver readers already return1110 for missing/invalid diagonal or
reciprocal-tensor material sidecars and1113 for a missing/invalid reference
subpixel configuration. For these cold errors the solver initializer destroys the
candidate and releases its child resources. Python now recognizes these two
unanimous statuses with sp_eh_init_recovery_protocol1, detaches, and releases its
failed lifetime reservation. Correct/stage the material file and retry in a new
work directory on the same caller communicator. Error messages identify the
material input instead of reporting a generic solver rejection.

No Fortran routine, input validation, coefficient arithmetic, file format,
communication or solver binary changed in this correction. The existing audited
runtime is reused by recorded SHA256. Unknown/mixed statuses or an absent/wrong
recovery protocol do not enable retry. Handled coefficient/Bloch-mode returns are addressed below; arbitrary allocation,
IO-device, MPI or process failure remain outside the bounded guarantee.

The existing E/H recovery tool accepts --material-input-mode diagonal|tensor|subpixel
and --material-input-admission. GNU MPI1/2 repeat missing, invalid-header, truncated
and nonfinite-file rejection twice each through the public Python constructor;
candidate/attachment resources are released and the parent communicator remains
usable. Valid material initialization then completes32 steps and direct restore.
Compare against the pre-change Python package with the same solver binary; this
is an input/lifecycle regression, not a new material-convergence certificate.

### Retrying handled coefficient/special-mode initialization errors

Python also recognizes existing solver cold returns1111 (nonfinite diagonal/
tensor LD update coefficients),1112 (nonzero Bloch phase with the unsupported
off-diagonal tensor mode), and1114 (invalid subpixel update coefficients).
Only unanimous codes with initialization recovery protocol1 release the failed
lifetime. The solver initializer already deallocates its candidate, releases its
send/receive cache before its owned MPI children, and returns; the correction is
Python's classification and diagnostic text, not a solver numerical change.

Correct the input/mode and use a new work directory on the same communicator.
This does not implement Bloch/off-diagonal tensor evolution, silently set Bloch
phase to zero, relax finite checks or change coefficient formulas. Files already
written during rejected initialization remain in the rejected work directory;
the rejection is not a promise of filesystem rollback. Unknown/mixed statuses,
missing protocol, trapped floating-point exceptions and allocation/device/MPI/
process failures remain outside this bounded retry guarantee.

--coefficient-admission in the existing E/H tool uses material-input-mode
(diagonal,tensor,subpixel) to exercise1111/1112/1114. Finite but oversized test
rates provoke the existing coefficient guards; a solver periodic test input
provokes the existing unsupported Bloch/tensor return. GNU MPI1/2 repeat each
failure and then run valid32-step/direct-resume controls. Reader rejections are
combined in the same fixture. The old SHA256-verified solver binary is reused;
no solver rebuild, extra numerical calculation in Python or new collective.

### Retrying handled early grid/boundary/time errors

FDTD updates of E and H already returns1101 for inconsistent cell/spacing/grid-count input,
1102 for incompatible periodic/boundary faces,1103 for an invalid duration/
step-count/time-step specification, and1104 when dt is not below its existing
CFL bound. Python now releases the failed initialization reservation for these
recognized unanimous returns only when initialization recovery protocol1 exists.
The solver cold candidate/owned-child cleanup remains the same. Correct the
input and retry with a new work directory on the same caller communicator.

No new input restriction, automatic repair, CFL threshold change, solver binary
change or collective is introduced. This aligns Python's error classification
with existing solver behavior. Unknown/mixed statuses, missing protocol and
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

FDTD updates of E and H initialization uses `allocate(...,stat=...)` for the material-ID and
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
allocated state; solver initialization cleanup releases it. One owning-field-
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

solver initialization detects STAT failures for three material vectors
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

### Memory for directional LD material initialization

For diagonal and off-diagonal LD permittivity tensors, Fortran checks `STAT`
when allocating arrays that retain material input. If one rank fails after
partial allocation, the field group agrees on1120 before broadcasting array
data. Partial arrays and owned communication resources are released. Python
returns to its pre-initialization state so the same input can be retried after
memory is released. Invalid input files still return1110.

The existing startup status broadcast is replaced by a field-group maximum
reduction; the number of collectives does not increase. Time-update arithmetic,
step communication and checkpoint formats are unchanged. Memory-budget excess
remains a warning; this error concerns actual required allocation failure.
It does not promise recovery from OS termination or arbitrary workspace failure.

The test-only `--material-allocation-failures` option injects failures after
partial diagonal, off-diagonal tensor and pole-coefficient allocation.
Production libraries contain no fault-injection hook.

### Allocation failures while preparing LD update coefficients

Fortran also checks allocation failures for the LD material list, six common
update-coefficient arrays and three directional coefficient arrays. Partial
failure skips coefficient evaluation and LD-state initialization that would
use those arrays. The field group agrees on1120, releases partial arrays and
owned communication resources, and returns Python to its pre-initialization
state. The same input can be retried after memory is released. Actual allocation
failure is an error; exceeding a planning memory budget still only warns.

Directional coefficients replace their existing status broadcast with a group
maximum reduction. Isotropic LD carries its status in the existing material-array
agreement. No startup or step collective is added. Coefficient formulas, arithmetic
order, time placement and checkpoint formats are unchanged. Nonfinite directional
coefficients retain1111, distinct from allocation failure1120.

The test-only `--ld-coefficient-allocation-failures` option fails one rank during
material-list, common2D/1D coefficient and directional-coefficient allocation.
Production has no injection hook. Tensor propagation workspaces, subpixel
reallocation, observation/image workspaces and OS/MPI failure recovery remain
outside this bounded support.

### Allocation failures for off-diagonal LD grid coefficients

The per-grid coupling weights, diagonal coefficients and material indices used
by off-diagonal LD also return1120 on allocation failure. Only allocation moves
earlier in initialization; value initialization and coefficient arithmetic keep
their original position and order. Allocation precedes the existing field-group
material-array agreement, so rank-local partial failure returns collectively
before use without adding communication. Partial arrays and owned communication
resources are released, preserving Python's pre-initialization state for retry
after memory is released. Budget excess remains a warning.

This applies to ordinary off-diagonal LD grid coefficients. Subpixel smoothing
allocates its own arrays and retains its separate behavior and coverage. Tensor
propagation workspaces, observation/image workspaces, OS termination and MPI
failure are outside this support. Formulas and checkpoint formats are unchanged.

The test-only `--tensor-grid-allocation-failures` option injects failures in the
5D coupling weights,4D diagonal coefficients and4D integer material indices.
Production libraries do not contain the fault-injection hook.

### Workspaces for off-diagonal LD time updates

Fortran allocates the large tensor-solver workspaces once at entry to an update
call and reuses them during iteration. The call owns this storage; physical
state, checkpoints and restore candidates do not contain it. Subpixel cases
without off-diagonal coupling do not allocate unused coupling buffers.

For Python stepwise calls, allocation status joins the existing field-group
input agreement. Partial failure returns24 on every rank before modifying E/H,
LD, clocks or outputs, releasing all partial scratch allocations. Python raises
`RuntimeFailure` while `field.failed` remains false. Release memory and retry
with the same J. This is distinct from `NumericalError`. No collective, formula,
update-order or checkpoint-format change is introduced. Budget excess still warns.

Coverage is limited to explicit tensor-update scratch allocation. It does not
promise recovery from physical-state backup allocation, compiler-generated
temporaries, arbitrary memory exhaustion, OS termination or MPI failure.
Standalone SALMON process execution stops on failed workspace allocation.
Observation/image allocation handling is not expanded indiscriminately.

The test-only `--tensor-workspace-failures` option fails one rank during drive,
solver and coupling-buffer allocation at early, intermediate and later updates,
checking unchanged state and retry. Production has no fault-injection hook.

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

This support applies to the SALMON batch FDTD path. The current FDTD updates of E and H
coupling API's unsupported-diagnostic return15 is unchanged. solver protected
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

This is the SALMON batch path;FDTD updates of E and H's unsupported-diagnostic return15 is
unchanged. Missing/partial spectra are not reconstructed or retried after
finalization. Check trajectory/energy/window coverage before analysis. This
optional-output policy does not silently accept failed required restart saves
or guarantee recovery from allocation/filesystem-wide failures.

### Required saved solver state-save failures

Unlike optional diagnostics,`sp_eh_checkpoint(...,0)` returns nonzero on save
failure. It uses `status='new'` without replacing existing files and agrees
failure within the owning group. Detached payloads leave live field/LD state
and clocks unchanged. After catching the error,the calling program can repair the destination
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

This is per-file staging. Binary payload publication waits until sidecar writing/close succeeds, but there is no atomic transaction across files/ranks or fsync durability guarantee. A failed phase can leave mixed generations; do not restart from that set. Retain a separately completed restart directory. This is distinct from FDTD called directly from Python protected restoration using independent candidates.

### Batch restart-directory preparation

An existing `data_for_restart` directory returns success status zero; symlinks to valid directories are also accepted. Directory-creation failure replaces the previous post-creation barrier with one owning-group status agreement, followed by nonzero termination on all participating ranks rather than a root-only stop. A conflicting existing file or symlink target is left intact. This handles creation of a required destination; it is not a writability preflight or a guarantee that later writes succeed. Writers check actual I/O failures. Ordinary-step communication and numerical updates are unchanged.

### Legacy batch EM restart-read errors

Text and real/complex binary readers require `io_status`. All owning-group ranks must call them in the same order, initializing the phase status to zero. Open/read/close errors retain the first failure; files open explicitly for reading. Shared reads distribute status before payload and never broadcast unread data after failure. Rank-local reads agree errors within the owning group. Remapping packs the previous three metadata broadcasts into one status-bearing header and checks workspace allocation and payload reads before transfer. No ordinary-step communication is added.

Legacy batch initialization exits nonzero immediately after a failed read, before consuming invalid values for setup or evolution. Helpers retain caller MPI and file ownership and support an explicit retry after repair and resetting status to zero. Failed output buffers/partially initialized models are not valid restored states. This differs from solver protected restoration using isolated candidates. There is no omitted-status legacy behavior or compatibility fallback.

`tools/test_em_restart_read.py BUILD_JSON NEW_OUTPUT --batch` runs GS-free missing/truncated/invalid-text/directory faults, repaired retries, and actual LD+FDTD midpoint resumes. The legacy observation-text copy strips trailing blanks from old rows; distinguish that formatting difference from numerical differences. It does not relax equality of physical values or binary restart payloads.

When input and output share a restart directory, new MPI-layout metadata is deferred until all required old payloads have been read successfully. A failed read leaves that input set intact. This is not an atomic multi-file transaction covering subsequent saves; retain a separate completed restart directory.

## Calculation-library records (for extensions and development)

NativeEHField accepts the existing in-process-tddft/1 calling program-library record and new
fast-hse-mpi/1 records declaring field_abi protocol1 with the same E/H timing
scheme. A new build_fast_hse_mpi build can therefore supply both FIELD_BUILD and
HSE_BUILD in the mixed_hse_models example; rank groups stay separate processes.
Its executable-build.json is the separate batch MPI FDTD executable record.
This does not add standard TDDFT/GS support to a fast library merely by changing
the record schema, and import does not compile or acquire a runtime.

</details>

## Performance mode for FDTD and embedded LD

After validating the input with the default mode, select the field mode explicitly at construction, identically on all field-group ranks. The choice is fixed for that object's lifetime.

```python
field = NativeEHField(
    context, calculation, build_record, work_directory,
    points=point_indices, grid=field_grid,
    execution_mode="performance",
)
A_next = field.advance(J_previous)  # physical charge-current density
```

Default `execution_mode="checked"` retains the existing path. Performance mode skips per-update full-array finite scans, rank-consistency comparisons of A/J, validation-only MPI agreements, automatic full-state copies and automatic rollback. It also applies to spatially owned `advance_owned`. Integration of A runs in Fortran.

Setup input/MPI/type/timing checks, build-time prohibited-solver audits and validated-library SHA checks remain. Minimal local lifetime/pointer-bound checks, physical halo/sample communication, tensor-solver convergence/breakdown handling and requested ordinary outputs remain. The caller owns valid current values, units, ownership and matching collective order.

Invalid values may propagate, and solver breakdown may terminate the whole job. There is no automatic state-preservation or in-place retry guarantee. Make explicit checkpoints for restart. Explicitly requested protected steps, energy queries and checkpoint operations retain their own validation/copy costs; they are not silently disabled.

This setting applies to the FDTD field and its internal LD update. It does not yet disable SBE/TDDFT/fast-HSE matter checks or outer coupling-driver checks; it is not an unchecked switch for the entire mixed calculation. The manifest records execution_mode.

Development fixture `tools/test_diagonal_field_energy_mpi.py` accepts `--execution-mode performance`. Omit `--observe` when measuring without energy diagnostics. Final checkpoint/resource reporting belongs to the fixture, not the performance step. Explicitly rebuild the field library containing `sp_eh_step_performance`.

## TDDFT performance mode

`NativeTDDFTBackend(..., execution_mode="performance")` selects unchecked
updates for the caller-managed MPI route. The default is
`"checked"`; the backend and returned state's mode are read-only. Supported
verification profiles are periodic fixed-ion scalar PZ, TB-mBJ, Libxc-PZ/PBE
with transverse supplied A/c and middlepoint propagation. Single-state, shared-model and heterogeneous-registry routes are supported.
Process-isolated workers, fast HSE and molecules are not covered by this option.

Initialization, library SHA256 authentication and explicit checkpoints remain
checked. Normal updates omit Python endpoint/finite/clock agreement, directory
admission votes, Fortran argument agreement, full-state backups and recoverable
numerical scans. Physical solver communication and update order remain. Shape,
capacity and lifetime bounds remain local checks. Caller ranks must issue valid
matching operations; nonfinite values can propagate and native failures can end
the MPI job. A returned update error marks the Python state failed; it cannot be
retried under the checked-mode state-preservation promise.

Explicit norm queries, observation, protected operations and output retain their
checks and costs. The outer coupling driver and FDTD have independent settings.
The sample `mpi_coupled_tddft` accepts `--matter-execution-mode performance` to
select this TDDFT route only; its field/outer checks are retained. The C entry
`sp_tddft_step_performance(Aold,Anew,J_matter)` returns SALMON matter current;
Python changes the sign once to charge current. Aold is ABI-only in this route.

Small synthetic GS tests compare checked/performance/previous library for all
four XC profiles on 1/2 internal ranks. Current and measured electron-number
arrays agree exactly, including explicit in-memory restoration. They establish
neither real-material convergence nor general performance gains.

## Per-state modes in shared and registered TDDFT

`NativeTDDFTBackend(..., shared_model=True, execution_mode="performance")`
and `registry_model=True` set the default for states created by normal sessions.
For direct model use, `model.new_state(point_id, execution_mode="checked")` or
`"performance"` overrides that default per state. Omitting the keyword inherits
the model default. `registry.open_model(..., execution_mode=...)` can set each
registered model's default. Each created state's choice is read-only; sibling
states may use different modes. Shared registry reuse preserves each binding's
selected default instead of inheriting a previous binding's mode.

Performance state-token updates still check local handle generation/lifetime,
capacity and ABI bounds. They omit token/endpoint/model-selection diagnostic
votes, numerical scans and automatic backups. Required model/input/workspace
activation and physical solver communication remain. Checked-state operations
keep the previous collective validation. Explicit snapshots/restoration remain
checked and do not change the receiving state's mode. No automatic rollback or
retry guarantee is implied for performance failures.

Synthetic regression covers interleaved checked/performance states, independent
controls, PZ/TB-mBJ/Libxc-PZ/PBE, heterogeneous PZ↔TB-mBJ and Libxc-PZ↔PBE, and
1/2 internal ranks. This is not a real-material convergence or speed benchmark.

## HSE performance mode

`NativeFastHSEBackend(..., execution_mode="performance")` selects the
caller-managed MPI route; the default is `"checked"`. The choice
is fixed on each state. Single-state, shared-model and heterogeneous-registry
routes are supported. Serial subprocesses use the separate `FastHSEResponseBackend` option described below.

The performance update omits endpoint agreement, Python finite scans, electronic
and exchange-state backups, and numerical-recovery diagnostic votes/scans.
Initialization, library authentication, local lifetime/capacity/shape bounds,
explicit checkpoints and observations remain checked. Solver breakdown checks,
FFT/Hartree/exchange reductions and learning-readiness communication remain.
The mathematical compatibility of A with the retained symmetry little group
also remains: status34 rejects that input before mutation and permits retry.
Other returned failures mark the Python state failed, without an automatic
rollback guarantee; unhandled native failures can end the job.

Learning endpoint collection stays before energy/output handling in both modes.
An owner-local scheduling flag separates this ordering from recovery diagnostics;
it is reset after each update and is not an execution-mode selector in shared
process globals. Taylor/PT-CN, full/reduced k and internal1/2-rank synthetic
comparisons preserve current, energy and learning history, including explicit
in-memory restoration. These are not physical convergence or speed benchmarks.

`mpi_mixed_native --hse-execution-mode performance` requires an explicit
`--fast-hse-mpi-build`. It selects only HSE updates;
the field and outer coupling retain their independently selected checks.

## Modes of shared and registered fast-HSE states

A backend's `execution_mode` supplies its model default. Use
`model.new_state(point_id, execution_mode="checked")` or `"performance"`
to override per state; omit the keyword to inherit the model default.
`registry.open_model(..., execution_mode=...)` selects a registered model's
default. Session registry reuse preserves each material binding's chosen mode.
Siblings can mix checked/performance updates with separate electronic, exchange,
learning, A/current and output history. There is no process-global mode selector.

The performance handle update keeps local token generation/lifetime bounds and
the symmetry/A mathematical condition, while omitting diagnostic token votes
and automatic recovery copies. Physical activation of each model's input and
workspaces, solver communication and learning update order remain. Checked
siblings retain their full validation and recovery path. Explicit memory or disk
checkpoint operations retain checks/costs; restoring does not change the
receiving state's selected mode. Performance update failures still have no
automatic rollback guarantee.

`shared_fast_hse --execution-mode performance` demonstrates point-a in the
selected mode and point-b in checked mode, including a matching point-a restart.
With `--second-initial`, the same example uses two registered input models.
The usual `mpi_mixed_native --hse-execution-mode performance` also permits
shared HSE points; field/outer-coupling checks have independent settings.

## Periodic TDDFT in a separate process: performance mode

`TDDFTInitialState.new_state(..., execution_mode="performance")` and
`crystal.respond(initial, drive, execution_mode="performance")` select the
periodic A-driven worker mode. `TDDFTResponseBackend(execution_mode=...)` passes
the choice to each worker used by a response session; outer session/field checks
retain their own settings. Default is `"checked"`, fixed per state. Molecules
and fast-HSE subprocess workers are not covered by this option.

The child must accept `INIT_PERFORMANCE` and reply with
`SALMONPY_TDDFT_PERFORMANCE_V1`; an old binary is rejected, without fallback.
Mode is passed only in that child's environment, never installed in the Python
host environment, and does not select embedded-host behavior. Setup, executable
SHA256, protocol shape/sequence, capacity and explicit output checks remain.
The step omits Python/worker endpoint-current finite checks and automatic
reconstruction/replay after a worker failure. A failed performance worker is
closed; do not continue it under the checked-mode state-preservation promise.
Physical propagation, update order, requested electron-number output and solver
breakdown handling remain. No global collective is added.

Accepted history is retained for explicit `snapshot()` / `restore()` and output.
This worker's explicit restore still starts a candidate and **recalculates saved
history**; it is not direct wavefunction restoration. Snapshot validation and
replay-current comparison remain, and the receiving state's mode is preserved
even when restoring a snapshot made in the other mode. Finite/nonfinite output
validation still applies when explicitly snapshotting or saving results. Normal
GS copying/integrity rules are unchanged; this mode does not reduce startup I/O.

Examples `tddft_response` and `mpi_tddft` accept `--execution-mode performance`.
Use a matching explicitly built or distributed worker; import does not build it.

## Outer MPI performance mode

Select `ResponseSession(..., execution_mode="performance")` or `CoupledStepper.from_multiscale(..., execution_mode="performance")` on every rank. The default is `"checked"`. Backend execution modes remain separate choices.

Array shape, ownership, step indices, and communication order remain checked. The mode omits per-step A/J finite scans and a redundant pre-update coupling agreement. Callback exception agreements and status-bearing A/J transport remain necessary to keep custom collective backends synchronized. `step_fields` still checks field contracts. Ordinary steps do not promise automatic rollback. Explicit `step_protected(targets)` retains backups, rollback, and its agreements in either mode.

Outer performance mode is also tested with actual E/H FDTD, built-in LD, and two TDDFT points: 3/6 ranks, 1/2 field partitions and 1/2 ranks per material group. All 96 steps of A/E/j match uninterrupted checked execution exactly after explicit protected rollback/retry in checked and performance modes. Existing synthetic initial orbitals test integration, not material convergence or speedup.

For HSE, select `NativeFastHSEBackend(..., execution_mode="performance")` independently of `ResponseSession(..., execution_mode="performance")`. Select `NativeEHField(..., execution_mode="performance")` separately for the field. Explicit `step_protected(targets)` retains internal backups of field, LD and HSE states for rollback/retry after a returned exception. It does not recover a terminated process or broken MPI. Reduced-k calculations still require the applied field to preserve the declared symmetry.

### A/E/H transport performance mode

Pass `ResponseSession(..., field_contracts=contracts, execution_mode="performance")` to `PacketCoupledStepper`. This omits duplicate finite scans of transported A/E/H values while retaining shapes, units, positions and step-index contracts. The pre-H error agreement remains because H acquisition is collective. Explicit `step_protected(targets)` also restores field-adapter indices and failure flags. Custom-material `StepTarget` callbacks must cover their internal state and side effects.

Molecular performance mode also supports coordinated direct disk save/resume of total-field FDTD, own-vacuum reference FDTD and the molecule. Separate MPI jobs verify restored Fortran step indices and prohibit propagation during restore. See [molecular coupling](DIPOLE_RADIATION.en.md) for conditions and diagnostic output keys.

## Choosing process-isolated HSE

`FastHSEResponseBackend(..., execution_mode="performance")` selects a fixed per-state mode; the default is `"checked"`. Keep this route when the host needs isolation from Fortran process state, working-directory ownership or MPI/main-thread constraints. It uses one serial CPU child per material point. Prefer `NativeFastHSEBackend` for internal MPI or reduced-k symmetry. The serial process route already requires full k points: reduced inputs cannot initialize, and performance mode does not extend that physical scope.

Before any update, the worker must accept `MODE_PERFORMANCE` and reply `SALMONPY_HSE_PERFORMANCE 1`. Unsupported workers are closed and rejected, with no silent fallback. Build explicitly with `tools/build_fast_hse_worker.py`; never compile on import.

Performance mode omits per-step finite scans, old-A comparison and optional Fortran diagnostic backups through the existing `rt_step(..., checked=...)` interface. Initialization, executable SHA, reply structure/index, shapes, capacity, lifetime and explicit save/restore checks remain. Physical update and learning/output order are unchanged. Failed workers are retired without automatic respawn or history replay; no new operation terminates independent siblings or the caller, but continuation of the failed state is not guaranteed.

Direct restore retains the receiving mode and existing caller-owned checkpoint authentication/A-restoration contract. `mpi_mixed_native --hse-execution-mode performance` selects either the MPI HSE backend or the serial worker. The outer coupling mode remains a separate choice.


## Many-point measurements

Use the manual `tools/measure_mixed_performance.py` driver for paired checked/
performance measurements with real E/H + two LD materials + independent SBE and
shared-model TDDFT states. It runs nine MPI ranks immediately; use allocated
compute resources. Both modes use explicit protected updates at every step.
See MEMORY guides for the required fixture, measurement meanings and limitations.
The observation stride does not thin the every-step macro-point current comparison.


For a controlled field-observation test use `measure_mixed_performance.py --study
observation`: fixed64 points/1024 steps,nine ranks,two repetitions of both modes
and strides1/16. It runs jobs immediately. All point histories remain every-step;
compare sparse E/H by common step number,not row index. MEMORY guides explain
exclusive per-rank phases and participant-summed Python MPI call counts. These
measurements neither remove protected backups nor alter the coupled algorithm.
