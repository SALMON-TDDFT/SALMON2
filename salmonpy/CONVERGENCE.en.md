# Grid and timestep diagnostics

`salmonpy.convergence` separates offline planning from execution. Python handles
inputs, trial dispatch and records. Finite-value checks, electron-number error
statistics and observable comparisons run in the bundled Fortran kernel. Import
and planning never compile, download or start a solver. An older native library
without this ABI requires an explicit wheel update or kernel rebuild.

`plan_grid_scan(cell_angstrom, radii=..., material_id=...)` starts at the smallest
selected physical PP radius/3 and refines to0.1 angstrom by default. Finer limits
are allowed. Each species record requires radius_angstrom, an explicit definition
and pp_sha256. Automatic radius extraction is not implemented: never substitute
the radial table extent or mesh step. Existing ceiling-grid resolution reports
actual spacings and removes duplicate effective grids. Keep a small fixed k mesh.

`run_grid_scan(plan, evaluate, atol=..., rtol=...)` accepts a callback returning
value (e.g. Hartree/atom), scf_converged (boolean) and optional metadata. Hold the
observable, units, geometry, XC, PP, bands and SCF conditions fixed. Unconverged
SCF warns and continues with provisional evidence. Adjacent-grid comparisons
require explicit tolerances; no passing pair yields None, never automatic success
at the finest grid. This stage establishes no k-point convergence.

```python
from salmonpy.convergence import plan_dt_scan, run_dt_scan
plan = plan_dt_scan(material_id="silicon", duration_fs=3,
                    impulse=[1e-4,0,0], expected_electrons=8)
report = run_dt_scan(plan, lambda dt: initial.new_state(dt=dt),
                     electron_number=backend.measured_electron_number,
                     atol=0.0, rtol=1e-5)
```

This is the connection contract, not a new method already attached to SALMONpy
TDDFT states. The backend must provide a Fortran/C measurement of occupied
wavefunction/density norm including weights and volume. Returning configured
nelec cannot detect integration errors. Built-in native TDDFT measurement wiring
and real-material validation are a subsequent stage. Each trial uses an independent
state at t=0,A/c=0, with initial current and step(A/c) returning three-component
charge current. States are closed explicitly.

The grid and GS remain fixed. dt starts at0.2 au and halves; the physical duration
is rounded up to a coarse step, then shared by every trial. The same weak jump
in A/c is applied at each endpoint; integrated E equals minus the impulse vector,
independent of dt. x/y/z and arbitrary vector polarizations are supported.

Primary accuracy evidence is the largest fractional electron loss. The default
relative target is **1/100000**. Gain and absolute errors are reported separately;
acceptance uses the maximum two-sided deviation. By default the initial measured
norm is the reference; expected_electrons instead checks against a known count.
Exceeding the target does not stop a trial. NaN/Inf, nonpositive norm or numerical
exceptions stop only that trial immediately; finer trials continue. Failed traces
are excluded from comparisons. Runtime/I/O failures and cleanup errors are distinct.
Invalid input/backend contracts propagate rather than masquerading as instability.

The report separates nonbreaking_dt_au from norm_acceptable_dt_au. These are
observations for the tested duration, not unconditional stability or physical
accuracy certificates. Optional current_atol/current_rtol compare adjacent finite
current traces at common physical times and report current_converged_dt_au.
Norm conservation alone cannot establish phase/current accuracy. No passing
criterion yields None. Check k-point convergence and production conditions separately.

Use `python -m salmonpy.examples.convergence_plan NEW_JSON` for a self-contained
offline planning example. Its explicit toy PP radius/hash are placeholders, not
validated silicon material data; it performs no electronic calculation.
