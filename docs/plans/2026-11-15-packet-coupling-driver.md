# Explicit packet coupling driver (first supported schedule)

Implement PacketCoupledStepper for existing NativeEHField/MultiscaleField and
packet-mode ResponseSession. One-to-one registered macro points only, exact native
samples, no automatic interpolation or extrapolation. Native available ports:
A/H offset+1/2, E offset0; output j offset+1/2. Multiscale A offset0, E lag1;
output j offset0. State clock remains local n*dt. Packet session accepts explicit
output j offset0 or1/2; field driver enforces the specific phase before any advance.
Legacy step(A)/CoupledStepper unchanged. Initial feedback remains explicitly zero.

Contracts must match field positions, atomic units, point sampling, t0/clockoffset0,
and supported offsets/lags. Require H ABI on all field ranks before advancing.
Driver: field.advance(previous j), read A/E (+collective H when required), packet
exchange, replace feedback. Failure invalidates response, no global rollback.
Return records contain separate field/port/current clocks. Test A-only numerical
parity with legacy, custom E/H closed loop, invalid clock/position/ABI refusal,
collective faults and old regressions. No isolated dipole auto-conversion, no
magnetic constitutive law, no new convergence study.
