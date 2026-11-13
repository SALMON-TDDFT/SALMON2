# Input parameter lifecycle

This is the SALMONpy public interface rule. LD implements it now; other backends
must adopt it when their Python interfaces are introduced.

| Lifetime | LD example | Allowed change |
| --- | --- | --- |
| Material/model | plasma frequency, pole strengths/resonances/damping | Construct a new material |
| Evolution state | dt, initial E, initial P/J, integration method | Construct a new state |
| Call | E at the next step | Supply each call |
| Query | Frequencies for permittivity, current P/J | Query without advancing state |

One parameter has one owner; `step` does not accept material settings or dt.
`respond(field, dt=...)` is a convenience operation that creates a fresh state
each call, not an override of an existing state's dt.

`material.config` and `state.config` expose frozen resolved records. Sequence
values are tuples, and `dataclasses.asdict(config)` produces an independent copy
for records/JSON serialization. The state configuration records initial settings;
`state.time`, `state.current` and `state.polarization` report evolving quantities.
Defaults (units, background permittivity, initial zero P/J and integration method)
are explicit in the records. They are fixed implementation choices, not options.

```python
from dataclasses import asdict
print(asdict(material.config))
state = material.new_state(dt=0.1)
print(asdict(state.config))
# state.dt = 0.2        # AttributeError; create another state instead
# state.step(E, dt=.2)  # TypeError; E is the only step parameter
```

Numeric scalars must be real, finite numeric scalars; strings, booleans, complex
numbers and length-one arrays are not accepted as scalar settings. Pole arrays
must be nonempty, have identical one-dimensional shapes and nonnegative finite
values. The initial E and each step's E must have shape (3,). The batch field must
have shape (N,3), N>=1. Inputs are copied so caller mutation cannot alter the model
or stored field. Public settings have no setters; private underscore attributes
are implementation details and must not be modified.

Validation failure, native-library failure and nonfinite native output do not
commit a step. P/J, previous E and step count are preserved. Step time is implicit:
the new E is at `state.time + state.dt`. No automatic reset or interpolation is
performed. The API cannot detect a caller supplying the wrong physical time.

For future backends, restart loading must check model, units, state schema and
time-evolution settings before accepting state. Configuration records alone are
not checkpoints: per-pole P/J and previous E are also required. No checkpoint
save/load API is implemented in this preview.

Output selection must not alter the physical update. An observable requiring
history or online accumulation must be enabled before that history begins;
requesting it late must not imply that missing past samples can be reconstructed.

SBE: eigenvalues, momentum/nonlocal matrices, weights, lattice, electron count,
correction order and nonlocal-correction flag belong to immutable SBEInitialState.
`new_state` fixes dt, time0 and initial A; only the next endpoint A is passed to
`step`. Density, previous A and current belong to each individual state. Failed
validation or nonfinite native output does not commit a step. Independent states
never share mutable density matrices.
