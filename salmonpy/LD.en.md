# Lorentz–Drude response

`salmonpy.LorentzDrude` provides isotropic passive-pole permittivity and physical polarization P/current J driven by E. Use SALMON atomic/Gaussian units: angular frequencies and damping in inverse atomic time, time in atomic time, and corresponding atomic E/P/J. Convert eV/fs explicitly; pole strengths are dimensionless.

```python
import numpy as np
from salmonpy import LorentzDrude
material = LorentzDrude(plasma_frequency=0.5, strength=[0.8, 0.2],
                       resonance=[0.0, 0.3], damping=[0.01, 0.02])
E = np.zeros((101, 3))
E[:, 0] = 1e-3
result = material.respond(E, dt=0.1)
with material.new_state(dt=0.1, initial_field=E[0]) as state:
    P, J = state.step(E[1])
    saved = state.snapshot()
    state.step(E[2])
    state.restore(saved)
    P, J = state.step(E[2])
```

Sample zero in `respond` is the initial P=J=0 at time zero. `step` accepts E at the next endpoint and uses stored previous E with a fixed dt and trapezoidal update. P/J are summed over poles, not FDTD auxiliary Q. Histories are independent. Check time-step convergence.

The continuous equations are `dP/dt=J` and `dJ/dt=f*wp**2*E/(4*pi)-gamma*J-w0**2*P`. `permittivity(omega)` uses exp(-i omega t) and returns `epsilon_infinity+sum(f*wp**2/(w0**2-omega**2-i*gamma*omega))`. Real-axis poles are rejected. Set instantaneous background permittivity with epsilon_infinity (default one); anisotropy, active media, fitting and automatic unit conversion are not included.

## Ownership and memory

Coefficients are computed once per state and retained natively with its per-pole history. Normal steps pass only three E components and receive three P/J components each. `current`/`polarization` return copies of accepted totals. Full pole histories cross the boundary only for snapshots/restoration. Nonfinite arithmetic rejects the entire call without changing accepted time/history, including batch calls.

Use a context manager or `close()` for direct states. `respond` and MPI response sessions close states automatically; GC also releases abandoned handles. Closed states reject stepping and snapshots/restoration; last time/totals remain readable. Do not operate on one state concurrently. Native registry calls are serialized.

There is no arbitrary memory or state-count limit. Dimensions and allocation failures are checked. Major native arrays occupy approximately `192*n_poles+48` bytes per state, excluding metadata, Python storage and input/output arrays. Many idle states retain memory, so close them when finished. Coefficients belong to states because dt can differ between states.

Material objects remain data-only and pickleable after use. For evolving states use [snapshot/restore](SNAPSHOTS.en.md) with the same model and settings rather than pickle. This does not add an LD plugin to the public C material vtable. Native builds remain explicit, never on import.

## Instantaneous background epsilon_infinity

```python
material = LorentzDrude(plasma_frequency=0.5, strength=[1.0],
    resonance=[0.3], damping=[0.02], epsilon_infinity=2.25)
response = material.respond(E, dt=0.1)
D = material.displacement(E, response.polarization)
u_background = material.background_energy(E, H)
```

The parameter is finite and positive, default 1. Inspect config.epsilon_infinity;
the legacy background_permittivity attribute is an alias. Material JSON accepts
the same optional parameter. P/J remain the dispersive-pole quantities.
displacement evaluates D=epsilon_infinity*E+4*pi*P. background_energy evaluates
(epsilon_infinity*E^2+H^2)/(8*pi), mu=1, excluding pole storage and dissipation.
Supply simultaneous, spatially aligned (...,3) arrays. These observables and
continuous permittivity are computed in Fortran. The background is not a
fictitious pole and belongs to the immutable model/checkpoint identity.

Native FDTD uses existing epsilon_em and includes the background in both
semi-implicit E-update coefficients. Pole update arithmetic/order is unchanged.
The conservative CFL includes vacuum/PML; backgrounds below one tighten it,
whereas a slower material does not relax the vacuum limit.

For non-unit backgrounds, Impulse remains the legacy dispersive-pole current
kick: propagation and raw pole/observation data continue. Legacy normalized
optical spectra assuming background one are omitted with a warning. Use
PlaneWave with an explicit reference/normalization for optical response.
The A-endpoint LD adapter also returns pole current only; the host field must
include the background separately, and this case warns. Do not feed only pole
current into dielectric_function and interpret the result as full epsilon:
that observable assumes total charge current. material.permittivity returns
the full continuous LD-model permittivity directly.

Local synthetic checks cover MPI1/2 slab/Fresnel agreement at background2.25,
the zero-dispersion dielectric limit, and direct native checkpoint continuation
with two distinct backgrounds. They do not certify real-material accuracy.
