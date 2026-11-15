# Common A, E and H input

The local field interface separates field data, required input ports and explicit
legacy adapters. Crystals still use A/c; LD uses E. A length-gauge model declares E.
A user model can require any declared subset of A, E and H. H is field intensity,
not B; neither a magnetic constitutive law nor an H-to-B conversion is inferred.
This interface does not by itself add isolated-system TDDFT or magnetic coupling
to any native solver. MPI ResponseSession still transports A only.

`FieldValue(value, time=..., position=..., unit=...)` owns an immutable real finite
3-vector and position. Time is always atomic time; position_unit defaults to bohr
and may be explicitly declared otherwise. Sampling defaults to point and supports
the existing TimeSpec sampling labels. `FieldPacket(A=None, E=None, H=None)` keeps
absence distinct from an explicit zero. Components can have different timestamps.

```python
from salmonpy.field_input import FieldValue, FieldPacket, FieldContract, FieldAdapter

# Existing crystal/SBE state, with declared A/c endpoint timing:
contract = FieldContract(state.timing, units={'A': 'A/c:au'},
                         position=(0., 0., 0.), gauge='velocity')
adapter = FieldAdapter(state, contract, legacy_quantity='A')
packet = FieldPacket(A=FieldValue([0., 0., 1.e-5], time=state.time+state.dt,
                                 position=(0., 0., 0.), unit='A/c:au'))
j = adapter.step_fields(packet)
```

For an E-driven state, supply a TimeSpec with InputPort('E'), unit E:au and
legacy_quantity='E'. Use gauge='length' for a declared length-gauge model, or
'phenomenological' for LD. The helper does not choose the gauge for a backend.
Velocity requires A/c, length requires E; a declaration is not verification of
backend physics. No A differentiation, H/B conversion or unit conversion occurs.
Legacy return values are preserved (LD returns P,J; SBE/TDDFT return current).
Legacy adapters require canonical atomic units and endpoint sampling with no
clock offset/input lag. Existing backend timing/dt/time, when exposed, must agree.
Native states can still be called directly with step(A) or step(E).

## Multiport user models

Build a FieldContract from a TimeSpec whose input ports are A/E/H, an exact
quantity-to-unit mapping, a position and gauge='custom' (or the appropriate gauge).
Pass a state implementing step_fields(packet) to FieldAdapter without
legacy_quantity. Unit labels must match exactly and have no built-in conversion.
An H unit such as SI:A/m is a caller declaration; a custom backend owns all physical
interpretation. The checker does not certify equations or units consistency across
quantities. The example records E/H separately, rather than adding unlike quantities.

At accepted index n, each required input time is
`t0 + (n + clock_offset + port.offset - port.lag)*dt`. Existing half-step offsets
can represent staggered H; time checks allow a small floating-point tolerance
(capped at dt*1e-3). Positions/position-unit and sampling must match exactly.
The interface represents fields already evaluated at one model position, not a raw
Yee-grid stencil: interpolate explicitly before handing off. It never moves values
between grid locations. At n=0 a port with uses_t0=False is not required.

Validation completes before stepping. Missing fields, bad shapes, complex/bool/
nonfinite values, wrong units/times/positions are rejected. Undeclared supplied
components are omitted from the packet delivered to a multiport backend. Selected
legacy vectors are writable copies; stored packet data remains immutable.

## State and failure ownership

FieldAdapter borrows the state and never closes it. step_index defaults to zero;
when attaching to an advanced state, explicitly supply its accepted index. Do not
advance/restore the backend separately while using the adapter; reconstruct it
with the correct index afterward. For both legacy and multiport backends, exposed timing/dt/time are checked before
and after stepping. Exposed time must equal t0+(accepted index+clock_offset)*dt.
Backends omitting these attributes rely on the caller declaration; physics is still
the backend responsibility.
Validation failure leaves the adapter usable without invoking the backend.
A backend exception invalidates the adapter; no backend rollback is invented.
This is not an MPI or whole-coupled-step transaction.

```sh
python -m salmonpy.examples.field_inputs
```

The installed-wheel example verifies native LD/direct equality and custom staggered
E/H delivery. Local native SBE/direct parity is covered by tests. Isolated native
TDDFT length-gauge connection and MPI E/H transport are follow-up implementation,
not capabilities established by these tests.
