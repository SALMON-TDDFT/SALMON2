# Molecular dipoles and FDTD coupling

## Single-molecule two-way coupling

Pass two already-created FDTD groups and one isolated TDDFT state to `MolecularFeedback`. Call collectively in matching order. Supply the molecular state only on `total_root`, a local `field` on every field rank, and `None` elsewhere. `reference_root` is the leader of the vacuum field driven only by this molecule, without illumination, other sources or material. Grid, boundaries, sample positions and current-deposition footprint must match the total field; setup checks this before advancing. The caller creates and closes MPI, fields and molecular states.

```python
from salmonpy.dipole_source import DipoleMap, MolecularFeedback

mapping = DipoleMap(('molecule',), (10,), [[1.]], [8.])
coupling = MolecularFeedback(
    world, field, state=state, mapping=mapping,
    total_root=0, reference_root=group_size,
    policy='environment_scattered_only', observer=observer,
)
first = coupling.start()  # one field interval with zero current
for _ in range(steps - 1):
    frame = coupling.step()
    if world.rank == 0:
        records.append(frame.record())
```

`world` is the caller's MPI communicator; `group_size` includes each field group and its bookkeeping ranks. `state` is `None` except on the reporting rank. `observer` is each field's `energy_observer(...)`, or `None` everywhere when not observing. Save `first` too if needed. `MolecularFrame` is returned only on the reporting rank; peers return `None`. It retains fields, current, dipole, electron number, freshness-aware internal energy and optional box energy observations. `record()` returns a detached serializable dictionary. There is no automatic file output.

The current entry supports one molecule/point with laboratory-axis dipoles and midpoint E driving. Initial selected E must be zero; otherwise increase source separation or delay. The field clock leads the electronic clock by dt. Subtracting the own-vacuum field removes vacuum radiation reaction too. Do not advance participating states outside the coupling while using it.

An inconsistent accepted endpoint or exception in the root-only `step(validate=callback)` rejects the interval. Catch `StepRecoveryError` on all ranks and retry only when `original_preserved=True`. There is no automatic retry. Validators must be read-only and perform no MPI operations. Electron restoration uses a direct live savepoint without replay. Worker loss or an indeterminate restore returns original_preserved=False and forbids retry.

The complete executable example is [molecular_feedback.py](examples/molecular_feedback.py). It explicitly reuses an existing GS and a validated E-driven worker; there is no implicit SCF, build or download.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 6 \
  python -m salmonpy.examples.molecular_feedback FIELD_BUILD E_WORKER RECEIPT EXISTING_GS NEW_OUTPUT \
  --field-grid 2 1 1 --ld
```

Its grid and duration are C2H2 wiring-test settings, not converged parameters for arbitrary molecules.

## Direct bundle restart: total field, own-vacuum field and molecule

Call `MolecularFeedback.save_checkpoint(directory, history=...)` collectively after
acceptance. Both FDTD groups must be at step n and the molecule at step n-1.
The field indices are read from Fortran. Do not save an unaccepted prediction.

```python
# First job; retain records only on total_root.
first = coupling.start()
records = [first.record()] if world.rank == total_root else None
while coupling.step_index < 30:
    frame = coupling.step()
    if world.rank == total_root:
        records.append(frame.record())
coupling.save_checkpoint('new-bundle', history=records)
# Close states and exit this MPI job.
```

```python
# A new job creates the same FDTD/molecular objects and coupling.
# Do not call start(); load the saved step indices directly.
saved = coupling.restore_checkpoint('new-bundle')
records = saved['history'] if world.rank == total_root else None
while coupling.step_index < steps:
    frame = coupling.step()
    if world.rank == total_root:
        records.append(frame.record())
```

`salmonpy.molecular-feedback-direct/1` stores both fields' E/H, LD, PML and
current buffers per field rank in Fortran, plus the molecular direct checkpoint.
It also retains sampled A/E, the saved E used by the next molecular update and
output streams. Optional caller-owned JSON `history` must be a list with one item
per accepted FDTD step, supplied only on total_root. Restore returns field_step,
molecular_step, per-owner verified_steps and history on total_root; peers get None.

Require the same runtimes, GS/input, dt, capacity, materials, boundaries, grid,
sampling and MPI decomposition. Verify hashes and configuration, stage both
fields and the molecule independently, verify pending Fortran indices and zero
molecular evolutions, then publish. Candidate-preparation failures preserve live
states and are reported collectively. Arbitrary OS/MPI failure or fatal allocation
failure during publication is outside this guarantee. A new directory and shared
filesystem are required. Normal updates gain no communication. Original states and candidates coexist
during restore; allow additional memory and output staging space. This entry needs
the field library's sp_eh_step_indices and revision5 molecular worker evidence.
Loading performs no build, SCF or elapsed propagation replay.

The complete examples/molecular_feedback.py supports
`--stop-after 30 --checkpoint NEW_CHECKPOINT`, then a separate MPI invocation
with `--resume NEW_CHECKPOINT`. Keep --steps (capacity), inputs and rank layout
identical; use a new output directory each time. --ld --field-grid 2 1 1 uses
MPI6. Measured molecular scope is PZ, fixed ions and one internal worker rank.


## One-way radiation sources

Use physical charge dipoles p in laboratory axes, in e·bohr. `DipoleMap` conserves
sum(p) when distributing polarization over explicit macro cells:

P_i = sum_m weights[i,m] p_m / cell_volumes[i]

j_i = sum_m weights[i,m] (p_m(new)-p_m(old)) / (dt cell_volumes[i])

Volumes are FDTD cell volumes in bohr³, never the vacuum-containing TDDFT box.
Weights are nonnegative, sum to one per molecule and are not silently normalized.
IDs are distinct molecule strings and nonnegative integer field-point IDs in the
specified order. Different molecules may overlap the same macro cells. Rotate
molecular dipoles into laboratory axes and apply actual molecule counts explicitly.
For a homogeneous ensemble, P=n*<p>; density and orientational averaging are caller
inputs, not inferred from a pseudopotential or simulation box. Units are atomic:
time in au, j in e/(bohr²·atomic_time). No SI conversion occurs.

```python
from salmonpy.dipole_source import DipoleMap
mapping = DipoleMap(('molecule-1',), (10,20), [[.25],[.75]], [8.,8.])
j = mapping.current([[0.,0.,0.]], [[0.,0.,.001]], dt=.003)
```

The map copies inputs and returns readonly arrays. Static dipoles produce no
current. Their electrostatic field is not initialized by this radiation source.
This is the electric-dipole approximation; resolved molecular charge distributions,
magnetic/multipole radiation and local-field/self-field feedback are not included.

## FDTD updates of E and H connection

`drive_native_dipoles(field, mapping, dipoles, dt=...)` borrows a fresh
`NativeEHField`. Call collectively on every field-group rank, using identical
inputs; other groups do not call this function. A precomputed trajectory has shape
(samples,molecules,3) at times n*dt from zero. The function returns readonly E/j
histories and their times, plus molecule/point IDs. E samples include t=0 and j
samples label interval averages at midpoint times. solver macro-cell to Yee-edge
mapping and FDTD numerical updates remain unchanged; do not deposit on edges twice.

The first current interval must be zero, following solver initialization. If an
existing molecular trace changes immediately, explicitly repeat its initial
sample. This delays the source by one dt, and is not hidden inside the API:

```python
import numpy as np
from salmonpy.dipole_source import drive_native_dipoles
# p: physical charge dipoles at n*field.dt, shape (samples,molecules,3).
# field: initialized NativeEHField with matching point IDs/cell volumes.
delayed = np.concatenate((p[:1], p), axis=0)
result = drive_native_dipoles(field, mapping, delayed, dt=field.dt)
```

The field must have capacity for len(delayed)-1 steps, matching point order and
volumes. Source dt must equal field dt; resampling is explicit. The complete
trajectory is validated and agreed before the first update. calculation implementation failure does
not roll back an entire trajectory. Caller closes the field. NativeEHField still
requires its explicit build record, MPI managed by the calling program context, interior vacuum source
neighborhood and plane-wave input. Use zero pulse amplitude to isolate molecular
radiation. An existing incident field can instead remain present; outputs then
contain both contributions. Boundary and CFL checks remain FDTD responsibilities.

This connection is one-way: the emitted E is not fed back to TDDFT. It does not
resolve radiation reaction or the distinction between self and local fields.
For the experimental isolated E-worker, use total charge dipole (last three reply
values), not current-sign conventions from periodic TDDFT. See FIELD_INPUT.en.md.

## Examples and verification

```sh
python -m salmonpy.examples.dipole_radiation
```

This lightweight wheel example checks current construction and conservation only.
In a source checkout, `tools/test_dipole_radiation_mpi.py` performs actual solver
FDTD tests with an explicitly supplied library called in the Python process record and a new output
directory, using two or three total ranks. An optional NPZ from
`tools/test_isolated_worker.py` supplies recorded molecular dipoles. See
`docs/dipole-radiation-validation-2026-11-15.md` for measured scope. The short C2H2
trace is a connection test, not a converged molecular emission spectrum.

## Refinement status

A fixed-position analytic point-dipole test now compares the complete retarded
E field and source work. At dt=.0015 au, refining dx from 2 to 1 to .75 bohr
reduced relative field error from 46.5% to 8.0% to 4.3%; source-work differences
from continuum radiated energy decreased from 15.8% to 2.5% to 1.1%.
This is not measured closed-surface flux or production-level convergence.
Halving dt alone did not monotonically improve continuum agreement. Source-cell
E is strongly grid-dependent, so raw total E feedback remains unvalidated.
Use tools/test_dipole_refinement_mpi.py and tools/summarize_dipole_refinement.py
in a source checkout; details and limitations are in
 docs/radiation-refinement-validation-2026-11-15.md. No self-field subtraction or
vacuum radiation-reaction model is silently enabled.

The sample also reads the public isolated_response output directly:

```sh
python -m salmonpy.examples.dipole_radiation --archive ISOLATED_OUTPUT --cell-volume 8000
```

Replace 8000 with the physical volume of each intended FDTD cell in bohr³.
This verifies units, uniform sample times, polarization-current construction and
conservation; it does not launch FDTD. Full solver propagation remains in the
MPI regression script described above.

## Explicit molecular drive selection

`select_molecular_field` selects already sampled electric `FieldValue` data.
It does not run a two-way solver. Always specify a policy:

- `provided_drive`: retain the caller-selected drive and record its physical
  `meaning` in provenance. Supplying a reference is an error.
- `environment_scattered_only`: subtract only this molecule's own-source vacuum
  reference from the total environment field. External illumination and other
  molecules' fields remain. This also removes vacuum radiation reaction; it
  does not represent full radiative dynamics or spontaneous emission.

Both samples must use E:au, bohr, point sampling, exactly equal times/positions.
Reference provenance requires `source_ids=[molecule_id]`, `medium='vacuum'`,
`incident=False`. Both provenance objects must declare matching `footprint_id`,
`sample_map_id`, `own_source_history_id`. These IDs describe matching deposition,
observation and own-source history; they are caller declarations, not proof that
separate Maxwell calculations used identical settings. Missing or inconsistent
references are rejected, never replaced with zero or an implicit total-field drive.

The returned `MolecularDrive.field` and `.metadata` are immutable; `.record()`
returns a detached JSON-compatible record containing input/reference/selected E,
position, time, policy, molecule ID and provenance. No interpolation or units
conversion occurs. The pure additive regression sample is:

```sh
python -m salmonpy.examples.molecular_field_selection
```

This selector does not advance fields. The centered oscillator schedule and work
checks are separate development tests below; isolated-TDDFT feedback is next.

### Development midpoint feedback (molecular TDDFT not attached)

Molecular work resumed on2026-11-22. The source tool
`tools/test_molecular_feedback_mpi.py` couples one test oscillator to two
independent fields, subtracting only its own vacuum-source field. Oscillator,
current, work and residual numerics are Fortran. This is a development substitute,
not a molecular model or a public two-way TDDFT driver.

After one explicit zero-current hold interval, capture both fields, predict
endpoint fields, update the oscillator using their midpoint, restore the fields
and advance with the corrected current. Recheck the actual endpoint update and
work before accepting. Trial observation files are rolled back too. The fixture
uses four ranks (one rank per field plus control ranks); it does not demonstrate
large spatial MPI scaling.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 2 python tools/test_molecular_feedback_mpi.py external BUILD NEW_BASE --dt .003 --steps 64
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python tools/test_molecular_feedback_mpi.py feedback BUILD NEW_RUN --baseline NEW_BASE --dt .003 --steps 64 --schedule centered
```

BUILD is an existing audited FDTD library build.json; no implicit build occurs.
Add `--environment ld` to both commands for the synthetic LD environment. Its
own-source reference contains neither LD nor illumination. Baseline dt, steps,
environment and library must match before updates. `--schedule held` selects an
explicit old-algorithm control; it is not a fallback.

At dt=.003 au the old held schedule had a temporal work mismatch about3.60e-11 au;
midpoint feedback reduced it to about2.67e-27 au for that fixture. The synthetic
LD test retains nonzero scattered feedback with consistent midpoint work.
Report normalized energy/work residuals as well as absolute errors. This checks
sampled E/dipole work, not full-volume Maxwell energy plus boundary flux,
overall coupled second-order convergence or arbitrary media/multiple molecules.
The chosen policy removes vacuum radiation reaction and includes no spontaneous
emission. The following development path connects the isolated E worker; C2H2 wiring is tested below, with physical convergence still separate.

An explicit `replay` control feeds accepted currents through ordinary FDTD without
trials/rollback. Vacuum and LD electric-field histories match exactly, checking
that the trial mechanism preserves the existing field update for the same input.

Dipole-to-polarization/current mapping and own-source subtraction also use the bundled Fortran kernel. Wheels include it; source users explicitly build the lightweight kernels. Mapping needs no SALMON executable, and import never compiles or communicates.


### Isolated E-worker connection (development)

The same fixture accepts `--isolated-executable`, `--isolated-validation`,
`--isolated-bundle` and `--isolated-input` to replace the oscillator with
`IsolatedTDDFTState`. Supply a validated isolated E worker and an existing GS
bundle. A periodic A worker is not sufficient. No SCF/download/build is implicit.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python tools/test_molecular_feedback_mpi.py feedback BUILD NEW_RUN --isolated-executable E_WORKER --isolated-validation VALIDATION_JSON --isolated-bundle EXISTING_GS --isolated-input rt.inp --dt .003 --steps 64
```

The first field interval holds zero current. Electronic time zero follows that
hold; output records `electronic_time` and the dt offset from field time. Initial
selected E must be zero to match the worker's E0; increase source separation or
delay otherwise. Fortran averages selected left/predicted-right E and passes that constant interval
field to the existing worker without changing its integration order. Fortran also
converts the dipole difference into charge current.

After field rollback and corrected-current advance, validate the actual selected
endpoint against the predictor. A mismatch rolls fields/electrons back and rejects
the interval. Native field rollback tokens are single-use: retire the predictor
token and capture another before the accepted attempt. Electronic cancellation
uses a direct live savepoint without replay.
`--inject-reject-step 30` rejects once and retries for development validation.
Other errors abort this MPI fixture; this is not a general recovery API.

In addition to the labelled synthetic worker, real fixed-ion C2H2 PZ TDDFT has
been connected. Electronic grid24×24×32, spacing0.35Å; one authorized GS calculation
is saved and reused. Maximum reference/E-worker dipole difference is about
3.86e-12 e·bohr. Four-rank vacuum/synthetic-LD tests (one field rank plus a control
rank per field) compare independent E-driven responses. The short64-step run at
dt=.003 au checks wiring, not converged optical response. Electron number, quantum work and inner-box six-face energy balance were subsequently verified as described below.

An MPI-linked single-process worker can hang if it inherits its host's MPI
connection. Child environments therefore omit OMPI/PMI/PMIX/PRTE data and rank
identities. Parent MPI and scheduler allocation IDs remain untouched. This startup
path is tested on this Mac/OpenMPI, not certified for all sites. The current public two-way entry is `MolecularFeedback`, described above.

E3 also records measured electron number and checks it during replay. Fortran diagnoses loss: the endpoint-control C2H2/LD64-step maximum loss fraction about9.95e-15 versus target1e-5. Retry histories including electron number match continuous execution exactly. This short-run diagnostic is not optical convergence evidence.


Internal-energy/work observation uses E3 and the existing energy cadence. Only
fresh endpoints participate; no propagation/eigen-energy frequency is changed.
Fortran separately accumulates held-endpoint work actually seen by TDDFT and
midpoint-field work. In the historical endpoint-control comparison with source amplitude .1,
the held-work relative residual is about0.897%, while midpoint work differs by
about28.49%. The default weak run is dominated by separately measured zero-field
energy drift; it cannot establish energy accuracy. No automatic background
subtraction, energy stop gate or full-coupled conservation claim is made.


### Molecular/field midpoint timing

The development isolated coupling defaults to midpoint driving.
Fortran averages selected left/predicted-right E after own-source subtraction.
Only accepted right E becomes the next left boundary; cancellation preserves both
quantum state and previous left E. Validate the midpoint from accepted fields.
No automatic endpoint fallback. The endpoint control remains a historical comparison record. Standalone IsolatedTDDFTState.step(E) and the
Fortran solver remain unchanged.

`drive_time` uses the field clock: interval center for midpoint, right boundary
for endpoint. Electronic state.time is the evolved interval endpoint. Record
work_applied separately from work_selected_field_midpoint. Total-field and
own-vacuum-reference work are also separate: selected-work alignment does not
prove total electromagnetic conservation.

For the C2H2 PZ amplitude.1 comparison, applied work equals selected-field midpoint
work numerically; internal-energy residual is about0.695%, versus28.49% with the
endpoint control. Do not subtract zero-field drift. Injected rejection/retry
matches every continuous history. The inner-box Maxwell/LD/six-face balance is described below.

## Energy balance through all six box faces (development)

`NativeEHField.energy_observer(lower, upper)` selects an inclusive box in zero-based grid indices, excluding sources, absorbing layers and the outer boundary. All field ranks call in matching order. Retain a pre-step backup, advance the field, then call `observer.sample_step()` before releasing the backup. This read does not advance the field. It uses one reduction within the owning field group, returning a dictionary on its leader and `None` elsewhere. Carry it in an existing exchange packet to avoid another application-wide communication.

Fortran computes six outward face fluxes, the field energy change, LD work and external-current work. LD work includes both stored polarization energy and dissipation. This does not estimate PML heating. Sideways scattered radiation is included in the six-face flux.

Verified profiles include vacuum, positive nondispersive background permittivity/permeability, scalar LD and diagonal LD, with nonnegative Ohmic conductivity. LD permeability remains one as required by the existing solver. Diagonal LD uses the x/y/z background components. Non-diagonal backgrounds use the separate constitutive exchange described below. Unsupported observation configurations warn and return `available=False`, while the calculation continues.

Fortran recovers the diagonal masses from the actual Yee update coefficients, not material IDs. With electric c1 and current cJ, c4=-4*pi*dt/cJ and epsilon=c4*(1+c1)/2; the implicit LD contribution cancels. Magnetic mass comes from the actual H curl coefficient. Interface ownership therefore agrees with the solver. The unit-background scalar-LD evaluation path retains the previous arithmetic.

At step n, the modified stored energy weights E(n)^2 by the electric diagonal and H(n-1/2)*H(n+1/2) by permeability. Separately report the conventional energy with centered H. No solver updates, checkpoint state or communication counts change.

The explicit development tool tools/test_diagonal_field_energy_mpi.py uses no GS and tests field MPI1/2, positive dielectric/magnetic backgrounds, scalar/diagonal LD and direct restart from step32. It executes calculations; it is not an input-only management-node checker. Off-diagonal background observation includes the constitutive boundary exchange from the actual interpolation. LD subpixel observation uses the actual smoothed background and dispersive response.

Yee E/H samples differ by half a step. `field_modified_*` uses the cross-product energy preserved by the discrete update; `field_physical_*` uses time-centered H in the conventional energy. Both residuals are retained. Closure of the discrete balance alone does not establish continuum energy accuracy.

With the molecule driven after own-vacuum-field subtraction, add the removed vacuum self-field work to the field change, outward flux, LD work and molecular internal-energy change. Do not subtract zero-field drift. The short C2H2 development test gave an accumulated discrete field residual of about 9.1e-15 Hartree and a joint residual of about 5.9e-8 Hartree, about 0.70% of molecular work. Normalization by the much larger incident field would hide this matter error. This is not a general physical-accuracy certification.

Real FDTD domain tests also passed with two partitions separately along x, y and z (outer MPI6: two field ranks plus one bookkeeping rank per context). Molecular saved arrays exactly equal the one-field-rank results; six-face integrated quantities agree within reduction-order roundoff. The development tool `tools/test_molecular_feedback_mpi.py` accepts `--field-grid NX NY NZ` and requires `2*(NX*NY*NZ+1)` total ranks. Rejection/retry with the same two-partition grid exactly matches continuous execution, including every saved array and six-face record. Electronic rollback now uses a direct live savepoint. This is separate from multi-node or large-decomposition evidence.

With diagonal background permittivity, off-diagonal LD pole strengths are supported: work uses the complete component currents, including cross-component response. Off-diagonal background permittivity is accounted for using the actual background interpolation and a separate constitutive boundary contribution.

For non-diagonal backgrounds, `constitutive_outward_energy` is the outward-positive constitutive boundary exchange. It is separate from `outward_energy`, the six-face Poynting flux. Old/new field energies include the cross contribution; both residuals include its boundary exchange. The locally assigned cross energy need not be positive by itself.

The current Fortran/C step observation has 19 components, ending with constitutive exchange and Ohmic work; the summary has 12 components. Explicitly rebuild older libraries from the current sources. Only non-diagonal background observation adds interpolation halo traffic and preparation-failure agreement within the owning field group. Ordinary time evolution, application-wide communication and checkpoint state are unchanged. Discrete balance is not a physical convergence certificate.

LD subpixel smoothing uses the same 19-component budget, with diagonal weights recovered from actual update coefficients and cross energy/exchange from the spatially varying background interpolation. Complete LD currents enter `ld_work`. The observer does not substitute original bulk coefficients.

This checks discrete balance of the reference-frequency smoothed model. It is not evidence of broadband material accuracy or curved-interface grid convergence. Cells that fall back to original coefficients are observed with those actual coefficients. Inspect applied/fallback sample counts and reference errors in `salmonpy-subpixel-summary.json`. SBE/TDDFT are not smoothed.

`ohmic_work` is the nonnegative energy dissipated in each accepted step: dt times cell volume times the sum of sigma E_mid squared, where E_mid averages old/new edge fields. Fortran recovers sigma from the actual electric/current update coefficients after removing the implicit LD contribution. Only subtraction roundoff is treated as zero.

Keep Ohmic dissipation separate from `ld_work` (work on dispersive polarization, not solely irreversible heat), external-current work, six-face Poynting flux and constitutive exchange. Both residuals and the summary include all terms. For zero-conductivity configurations the extra evaluation is skipped, existing observations are preserved and `ohmic_work=0`. There is no additional communication or change to propagation/checkpoint state. Explicitly rebuild the field runtime and lightweight kernel together.

## Molecular TDDFT performance mode

`initial.new_state(dt=dt, execution_mode="performance")` selects performance mode for an E-driven molecular worker. The default is `"checked"`. Performance mode skips per-step finite-value diagnostics and automatic recovery. Shape, capacity, response type, energy freshness, and explicit save/restore validation remain enabled. `step_fields` still validates its field contract. Automatic retry after failure is not guaranteed.

Snapshots do not select an execution mode: the receiving state retains its mode, including the candidate worker used for direct disk restore. Compare against checked mode before use. This option does not put the surrounding FDTD coupling into performance mode.
