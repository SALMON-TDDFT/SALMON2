# Field components (j→A) and coupling to matter

[日本語](FIELDS.md) · Topic: `macro-field` · Paired revision: `2026-11-14.2`

`salmonpy.maxwell` closes the light–matter loop from the field side. Matter
components (TDDFT, SBE, user models) take A and return a current (A→J). Field
components take the current and return A at the next time (j→A). `couple`
alternates the two:

```
j_0 = matter(A_0)  →  A_1 = field(j_0)  →  j_1 = matter(A_1)  →  A_2 = field(j_1)  → …
```

## Conventions

- Atomic units, laboratory frame. A is SALMON's `Ac = A/c`, and E = −dAc/dt.
- **Interfaces exchange the electric (charge) current density j.** Matter
  components keep SALMON's electronic matter current J internally, with j = −J.
  The `crystal` protocol returns j (SBE and TDDFT states included); a backend
  that declares J_matter output is converted by `couple` in exactly one place.
  When you drive a field component directly, pass j to `step`.
- Timing: `state.step(j_n)` takes the current at n·dt and returns A at (n+1)·dt (checked mechanically; see "Time-layout declarations").
  Matter advances one step with the new endpoint A and returns the current at
  that time. This is the order of SALMON's `time_evolution_step` / `calc_emfields`.
- As in SALMON, **the current at t=0 is not fed back**; the first `step` is
  determined by the external field alone.
- Settings are immutable. dt and the step capacity are fixed by
  `new_state(dt=..., steps=...)`; inspect `state.config` for resolved values.
- Invalid input or a nonfinite result leaves the state unchanged.

## Components

| Component | SALMON counterpart | Description |
|---|---|---|
| `Transverse(external)` | `trans_longi='tr'` | A = A_ext; the current does not act on the field |
| `Longitudinal(external)` | `trans_longi='lo'` | Bulk longitudinal geometry: A_ind(n+1) = 2A_ind(n) − A_ind(n−1) + 4π j(n) dt² |
| `ThinFilm(thickness, external, epsilon=(ε1, ε2))` | `trans_longi='2d'` | Normal incidence on a thin film between dielectrics; returns the transmitted field in the film |
| `Multiscale1D(cells, spacing, vacuum, incident, boundaries)` | multiscale, `fdtddim='1d'` | 1D Weyl-gauge propagation with matter points in cells 1..cells |

`external` and `incident` accept `AcosPulse`, a sequence of them,
`FieldSamples(kind='A')` (uniform fields only, on the n·dt grid), or a callable
that maps a time array to an (N,3) array.

### Multiscale1D details

- Cell i is at x = i·spacing; `vacuum=(left, right)` counts vacuum cells.
- The incident wave A_inc(t − x/c) enters through the left boundary. The initial
  field is built exactly as SALMON's `incident()`: the left region is filled at
  −3dt, −2dt and −dt, then one update produces the t=0 field.
- Boundaries: `'abc'` (first-order Mur), `'pec'`, `'periodic'`.
- x is the propagation axis; an incident wave with an x component is rejected.
- Stability requires c·dt ≤ spacing; c·dt = spacing is dispersion-free in vacuum.
- `state.grid_A` is A on the whole grid; `state.grid_E` is SALMON's centered E
  (located at (n−1)·dt).

## Time-layout declarations (`salmonpy.time/1`)

Each component declares, through its `timing` attribute, where its quantities
live in time (`salmonpy.timing.TimeSpec`). Before stepping, `couple()` checks
the declarations with `salmonpy.timing.check_loop` and rejects mismatches; it
never shifts by half a step or converts silently.

- Component: `dt`, `t0`, `clock_offset` (the n-th time is t0 + (n + clock_offset)·dt; 0 or ½)
- Input ports: `quantity`, `offset` (0 or ½), `sampling` (`point` / `interval_average` /
  `centered_difference`), `lag` (how many steps before output n+1 the latest input lies), `uses_t0`
- Output ports: `quantity`, `offset`, `sampling`, `at_t0` (`physical` or `zero`)
- Informative: `history`, `integrator`

Rules: (1) equal dt and t0; (2) connected ports agree in quantity, effective
offset and sampling (the only conversion is J_matter → j, a sign change);
(3) the loop's total lag is at least 1; (4) if a consumer uses t=0, the
producer's t=0 output is `physical`. Both declarations are kept in
`result.timing` and in the output archive's `timing` metadata.

| Component | Input | Output |
|---|---|---|
| Field components (all four) | j: offset 0, point, lag 1, t=0 unused | A: offset 0, point, physical |
| `SBEState`, `TDDFTState` | A: offset 0, point, lag 0 | j: offset 0, point, physical |
| Native E/H Yee field (`coupled.NativeEHField`) | j: offset ½, point, lag 1, t=0 unused | A: offset ½, point, physical |
| Matter in E/H coupling (response-service states) | A: offset 0, point, lag 0, on the electronic clock (+½) | j: offset 0 (+½) |
| User component without a declaration | `crystal` protocol default (as above) | as above |

A matter component without a declaration gets the crystal default (output j)
and a warning. `CoupledStepper` (E/H Yee coupling) gathers every point's
declaration over MPI, places it on the electronic clock (field clock + ½), and
checks it against `NativeEHField.timing`; the result is `stepper.timing`.

## Usage

```python
import numpy as np
from salmonpy.fields import AcosPulse
from salmonpy.maxwell import Longitudinal, Multiscale1D, SPEED_OF_LIGHT, couple

pulse = AcosPulse(1e-3, .2, 60., polarization=(0., 0., 1.))
result = couple(material, Longitudinal(pulse), dt=.02, steps=4000)
# result.A: (steps+1, 1, 3); result.current: charge-current density (steps+1, 1, 3)

layer = Multiscale1D(cells=3, spacing=SPEED_OF_LIGHT*.02, vacuum=(50, 50), incident=pulse)
result = couple(material, layer, dt=.02, steps=4000, output_directory='new-run')
```

`material` is an `AInitialState` (for example `SBEInitialState` or
`TDDFTInitialState`). A single one gives every point its own independent state;
pass a sequence with one entry per point to use different materials. `couple`
closes the matter states on success and on failure.

Example: `python -m salmonpy.examples.field_coupling NEW_OUTPUT_DIRECTORY`

## Caveats

- The coupling is explicit. The velocity-gauge current has a diamagnetic term
  proportional to A, which limits stability to roughly ω_p·dt < 2
  (ω_p² = 4πn_e).
- Not yet supported: coordinate rotation (SALMON multiscale `rmat_ms`), 1D media
  with a background permittivity, 2D/3D Weyl propagation, and restart
  checkpoints.
- E/H coupling on the Yee grid (`salmonpy.coupled`, [COUPLED.en.md](COUPLED.en.md))
  is a separate entry with a different time layout (E at n·dt, electrons at
  (n+½)·dt).

## Validation

`salmonpy/tests/test_maxwell.py` checks that:

- the longitudinal plasma-oscillation period matches the leapfrog dispersion
- thin-film transmission matches the analytic −iω/(−iω+γ)
- 1D vacuum propagation (c·dt = h) matches the analytic wave within 1×10⁻⁵
  relative, with no residual reflection
- a one-cell 1D layer agrees with the thin-film model
- `Transverse` coupling with the bundled SBE kernel is bit-identical to
  `crystal.respond`
- `Longitudinal` coupling is bit-identical to SALMON's recurrence

Direct comparison with SALMON itself: `tools/test_longitudinal_salmon.py`
(needs an external SALMON). The same synthetic Si initial state is propagated by
SALMON with `trans_longi='lo'` / `'2d'` / 1D multiscale (one point) and by the
SALMONpy TDDFT worker (`tr`) plus a field component. Over 400 steps, A and the
current agree within 3×10⁻¹⁵ relative, i.e. to SALMON's printed digits; in 1D
the whole-grid A agrees as well. See `docs/maxwell-field-validation-2026-11-14.md`.

Uniform fields need at least N+1 endpoint samples, from t=0 through N*dt, for N steps.
An incident wave requires an `abc` left boundary. Other left boundaries require `incident=None`.
