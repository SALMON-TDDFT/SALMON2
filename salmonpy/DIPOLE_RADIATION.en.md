# Molecular dipoles as one-way FDTD sources

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
delay otherwise. Pass predicted endpoint E to the existing worker without another
Python midpoint average or changes to its integration order. Fortran converts the
dipole difference into charge current.

After field rollback and corrected-current advance, validate the actual selected
endpoint against the predictor. A mismatch rolls fields/electrons back and rejects
the interval. Native field rollback tokens are single-use: retire the predictor
token and capture another before the accepted attempt. Electronic cancellation
uses **accepted-history replay**, not direct/no-recomputation restoration.
`--inject-reject-step 30` rejects once and retries for development validation.
Other errors abort this MPI fixture; this is not a general recovery API.

In addition to the labelled synthetic worker, real fixed-ion C2H2 PZ TDDFT has
been connected. Electronic grid24×24×32, spacing0.35Å; one authorized GS calculation
is saved and reused. Maximum reference/E-worker dipole difference is about
3.86e-12 e·bohr. Four-rank vacuum/synthetic-LD tests (one field rank plus a control
rank per field) compare independent E-driven responses. The short64-step run at
dt=.003 au checks wiring, not converged optical response. Electron number, quantum
work and full Maxwell energy conservation remain unvalidated.

An MPI-linked single-process worker can hang if it inherits its host's MPI
connection. Child environments therefore omit OMPI/PMI/PMIX/PRTE data and rank
identities. Parent MPI and scheduler allocation IDs remain untouched. This startup
path is tested on this Mac/OpenMPI, not certified for all sites. A general public
two-way TDDFT driver remains unavailable.
