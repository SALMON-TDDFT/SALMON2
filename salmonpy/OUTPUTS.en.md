# Reading and saving calculation outputs

[日本語](OUTPUTS.md)

`save_arrays` stores named arrays with explicit metadata in a new output directory.
`OutputArchive.open` verifies the archive; its `array`, `table`, `cube` and
`decomposed_current` methods read indexed data. These are postprocessing operations:
no solver is started, downloaded or compiled.

## Charge current, magnetization and spin-current dynamics

The existing Fortran solver writes `*_rt.data` and, for SOC, `*_rt_spin.data`.
`Jm_x/y/z` in the first file are matter current. Charge current is `j=-Jm`.
The Python `state.step(A)` and `state.current` already return charge current; do not negate them again.
Set `&analysis out_rt_spin_step=1 /` in the input for per-step spin observations.
The default stride is100; native rules also output steps0,1 and the final step.

```python
from salmonpy.outputs import read_spin_dynamics
spin = read_spin_dynamics('material_rt_spin.data')
t = spin.time
m = spin.magnetization       # (sample, x/y/z)
js = spin.spin_current       # (sample, spin polarization x/y/z, flow direction x/y/z)
js_zx = js[:, 2, 0]          # z-polarized spin flowing in x
print(spin.time_unit, spin.magnetization_unit, spin.current_unit)
```

Magnetization retains native `mag_x/y/z`, including component-sign conventions;
it is not converted to a magnetic moment in Bohr magnetons. The integrated density-matrix formulas are
`Mx=Re(rho12+rho21)/2`, `My=Re(-i*rho12+i*rho21)/2`, `Mz=Re(rho11-rho22)/2`.
The reader does not reinterpret these in another spin-coordinate convention. Spin current retains
the native Pauli-matrix definition, without an extra hbar/2 factor or sign conversion.
`pauli_current` retains all four0/x/y/z channels with axes `(sample, Pauli index, flow direction)`.
All arrays are read-only. Units come from the native headers without inference or conversion.
Physical observables are computed by Fortran; Python reads and organizes them.

## Save existing time series

```sh
python -m salmonpy.examples.outputs_roundtrip NEW_OUTPUT --rt-data material_rt.data --spin-data material_rt_spin.data
```

This extends the existing sample and does not run a calculation. It stores
`time_current`, `charge_current`, `time_spin`, `magnetization`, `spin_current` and
`pauli_current`. The sign conversion `j=-Jm` is explicit and uses NumPy's compiled
negate operation. Unit and definition metadata accompany the arrays.
Current and spin files can have different sampling times, strides and initial samples.
Their time axes remain separate: there is no row-number alignment or interpolation.
Restart output segments are not automatically stitched. Read completed output files,
not files being written concurrently by another process.

The Pauli identity channel matches native matter current Jm at matching times and units. Python charge current is j=-Jm.
See [NICKEL_SOC.en.md](NICKEL_SOC.en.md) for reading outputs and supported conditions.
