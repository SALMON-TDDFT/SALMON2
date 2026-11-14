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
