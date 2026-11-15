# MPI field packets (F4 transport subset)

Add optional per-point field_contracts to ResponseSession and collective
step_fields(mapping of transport point IDs to FieldPacket). Legacy step(A)
remains unchanged and cannot be mixed with the declared multiport session.
Contracts (ports, clocks, units, positions, gauge) agree before prepare callbacks.
Only field leader supplies packets; validate all required ports before broadcast
or any model update. Transfer plain value metadata and reconstruct immutable
packets; forward only declared ports. Backends explicitly implement step_fields
returning finite physical charge current density, using existing reduction/gather.
No implicit A/E conversion, interpolation, dipole normalization, or native E/H
coupling automation. Isolated dipoles require a separately explicit current model.
Test three ranks with different required ports, collective bad-time/missing-port
rejection before callbacks, callback failure invalidation, and legacy regression.
