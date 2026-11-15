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
j at the explicitly declared integer/half-step time from zero. Declare exactly one j output with point sampling and
offset zero or one-half; any backend output declaration must match. FieldAdapter also checks declared backend timing. Undeclared ports are omitted.
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

## Explicit time and space interpolation

`salmonpy.field_sampling.SampledField(quantity, times, positions, values, unit)`
retains immutable point samples with values[time,site,3].
`sample(time=..., position=..., weights=..., temporal='exact'|'linear')` returns
FieldValue. Supply nonnegative spatial weights summing to one and reproducing the
requested position. This validates a chosen stencil, not its physical suitability.
Exact sampling requires a stored timestamp; linear sampling requires bracketing
history. No extrapolation, nearest-time snapping, or unit conversion occurs.

`sample_packet(contract, step_index=..., histories=..., weights=..., temporal=...)`
uses each required port's own clock to construct FieldPacket. Histories and
weights are mappings keyed by required quantities. Only point sampling is
supported; interval averages are not relabeled. Required H needs an H history.

```sh
python -m salmonpy.examples.field_sampling
```

The example uses synthetic affine fields. A/E extraction is described below; the feedback driver supports only the schedules documented below. H extraction is described below.
Interpolation alone does not establish a self-consistent feedback scheme.

## Extract A/E samples from a coupled field

`NativeEHField.sample_fields()` and `MultiscaleField.sample_fields()` return
A/E SampledField objects after at least one successful step. This is a local
read without communication or updates. Samples are detached immutable copies;
rows follow field.point_ids and positions are in bohr. No H key is returned.

- Native E/H: A at (n+1/2)dt, E at n*dt; registered macro-cell centers.
- Multiscale1D: A at n*dt, E at (n-1)dt; material-cell x, y=z=0.

Retain samples over time and concatenate times/values with identical quantity,
units and positions into SampledField for explicit interpolation. Initial,
closed or failed fields are rejected. mpi_coupled_tddft saves field_samples.npz
with A/E, A_time/E_time, positions and point_ids, in atomic point-sample units.
This is not H data or a wavefunction restart. Native H extraction is described below; packet feedback is available through the limited schedules below.

## Explicit collective native H read

With an explicitly built H-enabled ABI, every field rank calls
`magnetic = field.sample_magnetic_field()`. After the first step it returns raw
H at (n+1/2)dt, spatially averaged to macro points, as SampledField. This differs
from native observation-file H, which is temporally centered at n*dt.
Reading does not advance the field or depend on observation stride. Libraries
without sp_eh_sample_h reject collectively. Unlike local sample_fields(), this
method must not be called only by the field leader.
The coupling example's --sample-h adds H/H_time to field_samples.npz. This does
not enable magnetic material physics, H-to-B conversion or packet feedback.

## Explicit packet feedback driver

Construct PacketCoupledStepper(response, field) and call step() on every context
rank. The response needs field_contracts; field is supplied only on field ranks.
Mapped material positions are one-to-one and must match exactly in atomic units.

| Field | Input port timing | Returned j |
|---|---|---|
| NativeEHField | A/H offset 1/2, E offset 0; all lag 0 | offset 1/2 |
| MultiscaleField | A offset 0 lag 0, E offset 0 lag 1; no H | offset 0 |

TimeSpec t0/clock_offset are zero; dt matches and sampling is point. State-local
time remains n*dt; physical port/output time includes its offset/lag. Declared
backend output must agree. H capability, unavailable future samples and position
mismatch are checked before advancing. Initial feedback is explicitly zero.
Order is field(previous j), sample, material, next feedback. A runtime failure
invalidates response without global rollback. Results retain field_time,
current_time, per-component SampledField values, and J_charge.

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_packet_coupling NEW_OUTPUT
mpiexec -n 3 python -m salmonpy.examples.mpi_packet_coupling NEW_NATIVE --field native --build BUILD_JSON --ports eh
```

Models in this example are synthetic. Dipoles are not current density, and no
magnetic material/self-field physics is inferred. Legacy CoupledStepper remains.
The new driver requires explicit packet backends; it does not silently adapt
existing TDDFT/SBE providers. General interpolation/subcycling remain follow-up.
