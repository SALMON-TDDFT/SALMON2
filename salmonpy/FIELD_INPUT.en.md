# Common A, E and H input

The local field interface separates field data, required input ports and explicit
legacy adapters. Crystals still use A/c; LD uses E. A length-gauge model declares E.
A user model can require any declared subset of A, E and H. H is field intensity,
not B; neither a magnetic constitutive law nor an H-to-B conversion is inferred.
This interface does not by itself add isolated-system TDDFT or magnetic coupling
to any native solver. MPI transport uses the explicit contracts and step_fields API below.

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
E/H delivery. Local native SBE/direct parity is covered by tests. Isolated native connection and MPI transport use separate APIs and validation described below.

## Experimental isolated native E connection (source checkout)

The bundled standard runtime now has a separate validation-only
`SALMONPY_TDDFT_WORKER=E1` protocol: `INIT_E dt capacity E0`, `STEP_E Ex Ey Ez`,
and replies `SALMONPY_TDDFT_E_V1 step delta_dipole[3] total_dipole[3]`.
Values use atomic units; dipoles use e·bohr. The electronic dipole change and total
charge dipole match native `rt.data` conventions. They are not current density.
Do not infer a dipole derivative or divide by the vacuum-cell volume implicitly.

This first validation covers isolated, fixed-ion, unpolarized PZ, scalar
pseudopotentials, one rank and midpoint propagation, with zero initial E. The
field passed at step n is the E used by native `calcVbox` at index n. The batch
reference prepares that sample from -(A[n+1]-A[n])/dt; no hidden A reconstruction
occurs in the worker. AETRS, external dipole sources and magnetic coupling are not
supported. Periodic TDDFT continues to use the original A-only protocol.

Developer checks: `tools/test_isolated_reference.py` prepares a small C2H2 batch
reference; `tools/test_isolated_worker.py` compares the new executable, interleaves
two independent child processes through FieldAdapter, and checks replayed E
history. `tools/isolated_worker_driver.py` is a validation tool, not an installed
public loader or an atomic snapshot/restore API. A backend error closes its child.
The tests do not establish optical convergence or full solver reentrancy.

Use an explicitly built, source-audited executable; existing release runtimes do
not thereby gain E support. The public isolated loader and in-memory replay are now documented in
ISOLATED_TDDFT.en.md, requiring a local binary-bound validation receipt.
Disk replay is documented in ISOLATED_TDDFT.en.md; MPI packets are described below. Prebuilt catalog distribution remains pending. See `docs/isolated-worker-validation-2026-11-15.md` in the source checkout.

## Explicit MPI A/E/H transport

Construct `ResponseSession(..., field_contracts={point_id: contract, ...})` on
all ranks, then call `session.step_fields({point_id: FieldPacket(...), ...})`.
Only the field leader supplies the mapping; others pass None. Keys are integer
transport point IDs from the plan. Required ports, clocks, units, positions and
gauge agree before state creation. All inputs are checked before updates;
receivers reconstruct immutable packets and validate them again. Backend states
implement step_fields and return three-component physical charge-current density
j at integer n*dt from t=0. Declare exactly one j output with point sampling and
zero offset; any backend output declaration must match. FieldAdapter also checks declared backend timing. Undeclared ports are omitted.
The legacy step(A) remains available without contracts; modes cannot be mixed.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_field_inputs
```

The sample uses synthetic current models. An isolated TDDFT dipole is not a
current-density return value. No A/E conversion, temporal/spatial interpolation,
or dipole normalization is inferred: supply explicitly sampled fields matching
FieldContract. Backend failure invalidates the session without whole-session
rollback. Reported field-value bytes exclude packet metadata and are not measured
network traffic. Native field/coupling drivers do not automatically use this API.
