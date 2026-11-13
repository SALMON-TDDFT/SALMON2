# SALMONpy: standalone LD preview

The first usable Python interface is an isotropic Lorentz–Drude material response.
It calls SALMON's extracted Fortran coefficient and time-update kernels; no SALMON
executable, MPI, ground-state calculation or FDTD grid is needed.

From the repository root, with Python 3.10+, NumPy, CMake and gfortran available:

```sh
python3 -m salmonpy.build_native
python3 -m unittest discover -s salmonpy/tests -v
python3 -m salmonpy.examples.ld_response
```

Set `FC` to the compiler executable if needed. The build is local and explicit;
importing the package never compiles code. The native library checksum is checked
on first use and its handle is cached. Restart Python (including notebook kernels)
after rebuilding the library. This is a preview with wheel packaging, not a PyPI release. Repaired wheel files
can be installed with `python -m pip install /path/to/salmonpy-....whl`. See
[PACKAGING.md](PACKAGING.md) for platform targets, CI and verified status.

```python
import numpy as np
from salmonpy import LorentzDrude

material = LorentzDrude(
    plasma_frequency=0.5,
    strength=[0.8, 0.2],
    resonance=[0.0, 0.3],
    damping=[0.01, 0.02],
)
E = np.zeros((101, 3))
E[:, 0] = 1e-3
result = material.respond(E, dt=0.1)
# result.time: (101,), result.polarization/current: (101, 3)

state = material.new_state(dt=0.1, initial_field=E[0])
P, J = state.step(E[1])  # fields at t=0.1; returns physical P and J
```

All values use SALMON atomic/Gaussian units: angular frequency, damping and plasma
frequency in inverse atomic time, time in atomic time, E/P/J in the corresponding
atomic units. Do not pass eV energies or femtoseconds without conversion. Pole
strengths are dimensionless. Each pole obeys

```
dP/dt = J
dJ/dt = f * wp**2 * E / (4*pi) - gamma*J - w0**2*P
```

`respond` starts with P=J=0 at the first supplied E sample. Sample zero is the
initial condition, not an update. `step` takes E at the next time, using the stored
previous E. Time step is fixed for each state. Every call to `respond` starts a
fresh state; separate `new_state` objects do not share histories. Observed P/J are
summed over poles. They are not the auxiliary FDTD quantity Q.

`permittivity(omega)` evaluates the continuous model with exp(-i omega t):
`1 + sum(f*wp**2 / (w0**2 - omega**2 - i*gamma*omega))`.
The instantaneous background is fixed at 1; anisotropy, active media, fitting,
and unit conversion are not included. An exact real-axis pole raises ValueError.
Finite-step time propagation approximates this continuous response; reduce dt to
check temporal convergence. A single state must not be shared between threads.

The step adapter reconstructs a small native context each call; it is intended
for prototyping, not high-throughput grid coupling. Batch calculations reuse one
native context. No Python FDTD or crystalline solver wrapper is provided yet.
Crystalline TDDFT/SBE interfaces remain A-driven only.

See [AI_USAGE.md](AI_USAGE.md) for extension rules and the legacy FDTD boundary.

Input ownership, immutable resolved settings and failure behavior are specified in
[INPUT_CONTRACT.md](INPUT_CONTRACT.md). Inspect `material.config` and `state.config`
for all effective settings. Python 3.10 or newer is required.
