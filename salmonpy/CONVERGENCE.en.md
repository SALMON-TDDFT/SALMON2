# Grid and timestep diagnostics

`salmonpy.convergence` separates offline planning from execution. Python handles
inputs, trial dispatch and records. Finite-value checks, electron-number error
statistics and observable comparisons run in the bundled Fortran calculation routine. Import
and planning never compile, download or start a solver. An older Fortran/C library
without this ABI requires an explicit wheel update or calculation routine rebuild.

`plan_grid_scan(cell_angstrom, radii=..., material_id=...)` starts at the smallest
selected physical PP radius/3, with0.1 angstrom as the default investigation limit.
It stops earlier once adjacent SCF-converged grids meet the declared tolerance. Finer limits
are allowed. Each species record requires radius_angstrom, an explicit definition
and pp_sha256. Automatic radius extraction is not implemented: never substitute
the radial table extent or mesh step. Existing ceiling-grid resolution reports
actual spacings and removes duplicate effective grids. Keep a small fixed k mesh.

Use `run_current_grid_scan(plan, evaluate)` for the primary grid decision.
For every grid, evaluate returns scf_converged and a dt_report from run_dt_scan.
Only a completed trial at norm_acceptable_dt_au is used. Compare adjacent grids
with the finer grid as reference. The default accepted error is **5% or less**:

    relative RMS = sqrt(sum((j_coarse-j_fine)**2)) / sqrt(sum(j_fine**2)) <=0.05

The sums cover common physical times and all three current components. This
relative L2/RMS metric avoids pointwise division at zero crossings. Different
dt values must have nested sample clocks; align by subsampling, never implicit
interpolation. A zero reference has no relative decision; only an explicit
positive absolute RMS tolerance enables an absolute-only decision.

Hold material, PP/XC, k mesh, excitation and physical duration fixed. A grid
without a norm-qualified dt is unresolved and excluded from comparison. SCF
nonconvergence warns and continues with provisional evidence. A confirmed
adjacent SCF-converged pair with error<=5% ends the scan before0.1 angstrom if
possible. Set stop_when_converged=False for the full plan. No passing pair yields
None. K-point convergence is separate. Scalar energy comparison via run_grid_scan
remains supplementary; it is not the recommended current-based grid criterion.

```python
from salmonpy.convergence import run_current_grid_scan
# evaluate_grid uses the Fortran/C backend and returns
# {"scf_converged": bool, "dt_report": run_dt_scan_result}.
report = run_current_grid_scan(spatial, evaluate_grid, rtol=0.05, atol=0.0)
```

```python
from salmonpy.convergence import plan_dt_scan, run_dt_scan
plan = plan_dt_scan(material_id="silicon", duration_fs=3,
                    impulse=[1e-4,0,0], expected_electrons=8)
report = run_dt_scan(plan, lambda dt: initial.new_state(dt=dt),
                     atol=0.0, rtol=1e-5)
```

New explicitly built runtimes expose electron_number on process TDDFT states
and ordinary native/shared/heterogeneous-registry states. It is measured in Fortran
from occupied, k-weighted spin density and grid volume, not configured nelec.
run_dt_scan uses this property by default. A custom Fortran/C calculation implementation may instead
supply electron_number=callback. solver queries are collective: all internal MPI
ranks call in the same order. Queries are read-only, state-handle-specific and
measure restored density after restart. Supported configurations are crystalline
A-driven PZ/TB-mBJ/Libxc-PZ/PBE. Fast PR HSE and isolated E-drive measurement
are not connected yet. Legacy runtimes still support ordinary calculations;
missing diagnostics fail explicitly rather than substituting an input count.
Updating the wheel alone does not update an existing solver binary.

```sh
mpiexec -n 2 python -m salmonpy.examples.tddft_convergence BUILD_JSON NEW_OUTPUT --grid 16
```

The example retains GS in memory and explicitly disables GS restart output.
Each dt creates independent material state from the same GS. SCF nonconvergence
warns and continues with provisional evidence. This small low-k example does not
certify optical convergence. The solver's existing order-of-magnitude electron
count guard returns a recoverable numerical error in the recovery-capable API,
restoring the accepted state and ending only that trial. This is distinct from
exceeding the1e-5 diagnostic target.

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
Invalid input/calculation implementation contracts propagate rather than masquerading as instability.

The report separates nonbreaking_dt_au from norm_acceptable_dt_au. These are
observations for the tested duration, not unconditional stability or physical
accuracy certificates. Optional current_atol/current_rtol compare adjacent finite
current traces at common physical times and report current_converged_dt_au.
Norm conservation alone cannot establish phase/current accuracy. No passing
criterion yields None. Check k-point convergence and production conditions separately.

Use `python -m salmonpy.examples.convergence_plan NEW_JSON` for a self-contained
offline planning example. Its explicit toy PP radius/hash are placeholders, not
validated silicon material data; it performs no electronic calculation.

The5% relative-RMS current default is a grid-screening criterion, not a dielectric
or high-harmonic spectrum accuracy guarantee. Choose stricter `rtol`/`atol` and
check the target spectrum separately where required; the report records the
actual tolerances. The default numerical threshold is unchanged.
