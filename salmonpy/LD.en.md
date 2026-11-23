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

The continuous equations are `dP/dt=J` and `dJ/dt=f*wp**2*E/(4*pi)-gamma*J-w0**2*P`. `permittivity(omega)` uses exp(-i omega t) and returns `epsilon_infinity+sum(f*wp**2/(w0**2-omega**2-i*gamma*omega))`. Real-axis poles are rejected. Set instantaneous background permittivity with epsilon_infinity (default one); this scalar model is isotropic and passive, and does not convert units automatically. Use `TensorLorentzDrude` for direction-dependent response. Constrained fitting is a separate explicit entry described below.

## Ownership and memory

Coefficients are computed once per state and retained natively with its per-pole history. Normal steps pass only three E components and receive three P/J components each. `current`/`polarization` return copies of accepted totals. Full pole histories cross the boundary only for snapshots/restoration. Nonfinite arithmetic rejects the entire call without changing accepted time/history, including batch calls.

Use a context manager or `close()` for direct states. `respond` and MPI response sessions close states automatically; GC also releases abandoned handles. Closed states reject stepping and snapshots/restoration; last time/totals remain readable. Do not operate on one state concurrently. solver registry calls are serialized.

There is no arbitrary memory or state-count limit. Dimensions and allocation failures are checked. Major solver arrays occupy approximately `192*n_poles+48` bytes per state, excluding metadata, Python storage and input/output arrays. Many idle states retain memory, so close them when finished. Coefficients belong to states because dt can differ between states.

Material objects remain data-only and pickleable after use. For evolving states use [snapshot/restore](SNAPSHOTS.en.md) with the same model and settings rather than pickle. This does not add an LD plugin to the public C material vtable. solver builds remain explicit, never on import.

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

FDTD called directly from Python uses existing epsilon_em and includes the background in both
semi-implicit E-update coefficients. Pole update arithmetic/order is unchanged.
The conservative CFL includes vacuum/PML; backgrounds below one tighten it,
whereas a slower material does not relax the vacuum limit.

For non-unit backgrounds, Impulse remains the legacy dispersive-pole current
kick: propagation and raw pole/observation data continue. Legacy normalized
optical spectra assuming background one are omitted with a warning. Use
PlaneWave with an explicit reference/normalization for optical response.
The A-endpoint LD adapter also returns pole current only; the calling program field must
include the background separately, and this case warns. Do not feed only pole
current into dielectric_function and interpret the result as full epsilon:
that observable assumes total charge current. material.permittivity returns
the full continuous LD-model permittivity directly.

Local synthetic checks cover MPI1/2 slab/Fresnel agreement at background2.25,
the zero-dispersion dielectric limit, and direct saved solver state continuation
with two distinct backgrounds. They do not certify real-material accuracy.

## Constrained fitting and provenance

```python
from salmonpy.ld import fit_lorentz_drude, LorentzDrude
config, report = fit_lorentz_drude(omega, epsilon, poles=2, drude=True,
    epsilon_infinity=True, units="atomic_gaussian", optimizer="scipy", rtol=0.05,
    source_file="my_optical_data.txt", reference="Caller-supplied citation")
if config is not None:
    material = LorentzDrude.from_config(config)
print(report.as_dict())
```

Supply increasing nonnegative angular frequencies (strictly positive for Drude fits) and complex epsilon, or nk=(N,2)
columns with n,k. The n/k conversion assumes a nonmagnetic medium (mu=1); supply epsilon explicitly for magnetic media. The absorption convention is n+i*k with exp(-i omega t);
conversion to epsilon is solver. Units must be atomic_gaussian, without implicit
conversion. A dielectric_function result may be supplied directly as the first
argument; masked bins are excluded and the optical configuration is retained.
Select a polarization/component explicitly; no tensor averaging or data download
is performed. Experimental tables must be supplied by the user.

poles counts the Drude pole. drude=True fixes its first resonance to zero.
epsilon_infinity=True fits a positive background, False fixes one, and a positive
real fixes that value. initial accepts an LDModelConfig or LorentzDrude seed.
Returned coefficients use plasma_frequency=1, with oscillator weights in strength.
solver projected LM keeps strengths/resonances/damping nonnegative and the fitted
background at least1e-9. The objective is sum|delta epsilon/max(1,abs(epsilon))|^2.
relative_rms_error is the unweighted complex residual norm divided by the data
norm; rtol (default5%) classifies that quantity, not universal physical accuracy.

starts (default4) selects deterministic initial guesses; max_iterations (default300)
is per start. This is local optimization: compare caller-selected pole counts,
seeds and bands to assess ambiguity and overfitting. An inaccurate or unfinished
optimization warns and may return a passive candidate; config=None means no finite
candidate. accepted and optimizer_converged are separate. No existing simulation
is stopped and no additional solver job is launched.

FitReport records band, both residuals, model/iteration settings, initial model,
software version, calculation routine elapsed time, citation, reference-file SHA256 and separate
normalized-array SHA256. source_file is a caller-supplied reference, not a claim
that this API parsed it. as_dict() returns independent JSON-ready data.
recommended_dt=0.1/fastest active pole/plasma/damping scale is an accuracy heuristic,
not a stability bound. Stability depends on the field grid/CFL and remains unknown.

The calibrated band follows config.valid_frequency_range through JSON/providers
and FDTD. Out-of-band carriers or broadband impulses warn; an in-band carrier does
not certify finite-pulse tails. material.check_excitation(spectrum_of_E, rtol=.01)
reports out-of-band raw FFT-bin squared-amplitude fraction, not absorbed energy.
It never blocks a running calculation. Run
`python -m salmonpy.examples.ld_fit --poles 1 2 3` for explicit pole-count comparison
on synthetic data, without real-material or global-optimum claims.

See [OPTICAL_DATA.en.md](OPTICAL_DATA.en.md) for the offline CC0 measurement and stored-candidate library.

### Optional SciPy optimizer

Install `salmonpy[fit]` explicitly and select `optimizer="scipy"`. SciPy uses bounded
`least_squares(method="trf", x_scale="jac")`; Fortran evaluates the complex residual,
initial guesses and final error/timestep diagnostics. SciPy controls optimization
and finite-difference Jacobians. No Python implementation of the LD physics is added.
The default is `optimizer="native"`; missing SciPy raises an
actionable error, without silent fallback or installation. Stored models need no SciPy.

For SciPy, max_iterations limits function evaluations per start (`max_nfev`); numerical
Jacobian calls are additional. The report records SciPy version, all start outcomes,
selected start, tolerances, optimality and termination reason. Local termination
does not guarantee a global optimum or acceptable physical accuracy.

## Common isotropic tensor configuration

`model.tensor_config` / `config.tensor_config` returns LDTensorModelConfig with one independent
component, epsilon_ij=epsilon*delta_ij. `LorentzDrude.from_config(tensor)`, FDTD materials
and LDMaterialProvider accept it. This applies across modes, not only3D E/H.
`permittivity_tensor(omega)` returns immutable complex(...,3,3) arrays; scalar APIs remain.
FDTD.tensor_materials and its run record expose the representation. FDTD called directly from Python now supports isotropic, diagonal and reciprocal off-diagonal tensors; standalone/MPI tensor material response is now implemented below.
See OPTICAL_DATA.en.md for oxide data, source assumptions and limits.

## Reciprocal diagonal/off-diagonal LD response

TensorLorentzDrude accepts a real symmetric3x3 strength matrix per pole, a symmetric
positive-definite epsilon_infinity background, and nonnegative scalar resonance/damping
per pole. Strength matrices are positive semidefinite; individual off-diagonal
entries may be negative. epsilon=epsilon_inf+sum(wp**2*F_q/(w_q**2-omega**2-i*gamma_q*omega)).
Fortran owns matrix validation, frequency response, P/J integration and background D/energy.

```python
from salmonpy.ld import TensorLorentzDrude, LorentzDrude

model = TensorLorentzDrude(plasma_frequency=1.,
    strength=[[[.12,-.04,.01],[-.04,.20,.02],[.01,.02,.30]]],
    resonance=[.4], damping=[.03],
    epsilon_infinity=[[2.,0.,0.],[0.,3.,0.],[0.,0.,4.]])
epsilon = model.permittivity_tensor([.1,.2])  # (2,3,3), lab xyz
with model.new_state(dt=.01) as state:
    state.step([.001,0.,0.])
    saved = state.snapshot()
    state.step([0.,.001,0.])
    state.restore(saved)
copy = LorentzDrude.from_config(model.config)
# LDMaterialProvider(model.config) accepts the configuration for MPI groups.
```

permittivity() does not silently reduce a tensor to a scalar. Select components explicitly.
Owned pole P/J and accepted E support direct snapshot/restore without replay. Rejected
nonfinite steps retain accepted state. Concurrent use of one state is not supported.
A-endpoint adapters return dispersive-pole current only; hosts own the instantaneous
background tensor. This is reciprocal passive response, not antisymmetric magneto-optics.
Automatic axis rotation remains unimplemented. Tensor fitting is described below.
Standalone/MPI material response and FDTD called directly from Python support reciprocal off-diagonal tensors without dropping components.
The FDTD updates of E and H implementation uses paired Yee edge/node means and their transpose; geometry is not restricted to a 3D-only material API.
Run `python -m salmonpy.examples.ld_response --anisotropic` for a synthetic example.


Three-drive TDDFT data and reciprocal tensor fitting

`analysis.dielectric_tensor` accepts three independent, simultaneous charge-current
and total-electric-field histories of shape `(3,N,3)`. Each frequency solves
`sigma=J E^-1`, then `epsilon=I+4*pi*i*sigma/omega`. SALMON matter current must
be explicitly negated before this call. No automatic A-to-E conversion or unit
conversion occurs. DC, weak excitation and rank-deficient runs are masked, never
filled with zeros. Inspect linearity, duration/window and convergence separately.

`fit_lorentz_drude_tensor` requires optional SciPy and an explicit initial tensor
model. It fits common pole frequencies/damping and PSD strengths `F=L L^T`, with
its initial SPD background fixed. Fortran computes matrices and residuals; compiled
LAPACK factors the initial matrices and SciPy optimizes. All nine measured components
remain in the objective, including antisymmetric residuals. This reciprocal passive
model cannot reproduce general magneto-optics. Accuracy and optimizer termination
are distinct report fields. An inaccurate candidate is returned with a warning.
A zero initial pole factor can be stationary: provide nonzero initial strengths.

The new extraction/fit path has synthetic numerical verification, not a newly
converged real-material TDDFT validation. Reciprocal off-diagonal FDTD called directly from Python is connected as described in FDTD_MATERIALS.en.md. A fitted candidate still needs accuracy, bandwidth, stability and grid checks.

```python
from salmonpy.analysis import dielectric_tensor
from salmonpy.ld import TensorLorentzDrude, fit_lorentz_drude_tensor

# time: (N,), charge_current and total_E: (3,N,3).
# First axis selects three independent impulse runs; last axis is lab xyz.
epsilon = dielectric_tensor(time, charge_current, total_E,
    current_convention="charge", field_kind="total",
    units="atomic_gaussian", window="none",
    metadata={"runs": ["impulse-x", "impulse-y", "impulse-z"]})
initial = TensorLorentzDrude(plasma_frequency=1.,
    strength=[[[.3,.02,0.],[.02,.4,0.],[0.,0.,.5]]],
    resonance=[.7], damping=[.08],
    epsilon_infinity=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
model, report = fit_lorentz_drude_tensor(epsilon, initial=initial,
    units="atomic_gaussian", rtol=.05)
# Inspect report before use; this does not start an FDTD calculation.
```

## Standalone performance mode

Pass `execution_mode="performance"` to `new_state(...)` or the batch response
entry (`material.respond(...)` for LD, `crystal.respond(...)` for SBE). The
default is `"checked"`; the choice belongs to each state and cannot be changed
during stepping. LD supports scalar and reciprocal tensor materials; crystal
SBE remains driven by A/c and returns charge current density.

Performance mode uses the same Fortran update arithmetic. LD omits per-step
pole backups and finite-value scans; SBE updates its density matrix directly
without a trial copy or per-step finite-value scans. Automatic rollback and
retry after a failed update are not guaranteed in this mode. Initial input,
library integrity and shape checks remain. Explicit `snapshot()` / `restore()`
and saved output retain their checks. Restoring physical state does not change
the execution mode. This option does not disable checks in TDDFT or external
coupling drivers.

The installed examples accept `--execution-mode performance`. Compare a checked
run before using this mode for a new configuration.


### Explicit damping floor for scalar fits

`fit_lorentz_drude(..., minimum_damping=0.002)` constrains every pole, including
Drude, in atomic units. Both the Fortran optimizer and optional SciPy optimizer
apply the same lower bound. Initial widths below it are projected onto the
allowed range; the input initial model is not mutated. The default is0, which
keeps the previous unconstrained calculation. The returned report records
`minimum_damping_au` and the initial model for reproduction. This is a model
assumption; a small residual on discrete samples does not validate narrow peaks
between them. It does not alter existing LD time evolution or loaded models.
