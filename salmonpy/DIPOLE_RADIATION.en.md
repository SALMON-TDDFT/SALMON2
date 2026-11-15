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

## Native E/H connection

`drive_native_dipoles(field, mapping, dipoles, dt=...)` borrows a fresh
`NativeEHField`. Call collectively on every field-group rank, using identical
inputs; other groups do not call this function. A precomputed trajectory has shape
(samples,molecules,3) at times n*dt from zero. The function returns readonly E/j
histories and their times, plus molecule/point IDs. E samples include t=0 and j
samples label interval averages at midpoint times. Native macro-cell to Yee-edge
mapping and FDTD numerical updates remain unchanged; do not deposit on edges twice.

The first current interval must be zero, following native initialization. If an
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
trajectory is validated and agreed before the first update. Backend failure does
not roll back an entire trajectory. Caller closes the field. NativeEHField still
requires its explicit build record, host MPI context, interior vacuum source
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
In a source checkout, `tools/test_dipole_radiation_mpi.py` performs actual native
FDTD tests with an explicitly supplied embedded library record and a new output
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
conservation; it does not launch FDTD. Full native propagation remains in the
MPI regression script described above.
