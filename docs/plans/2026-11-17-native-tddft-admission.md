# Recoverable embedded TDDFT admission

First bounded part of the approved native error contract: native (not only Python)
collective validation before initialization latches or accepted numerical state
are changed. Introduce a queryable rejected-before-update protocol. Status10:
invalid values;11: rank disagreement;12: unavailable lifecycle;13: capacity;
14: pending restore candidate. Every member of the bound solver communicator must
participate; process/MPI failures are excluded. Rejected initialization remains
retryable if native initialization has not started. A successful initialization
still owns the one-lifetime latch. Step checks reject before any RT mutation;
valid-call arithmetic/order stays unchanged. Kernel errors/STOP, OOM, I/O rollback,
whole mixed-step atomicity and HSE are not covered by this first subset.

Verify actual two-rank raw ABI bad calls, snapshot byte equality before/after
rejections, valid retry vs uninterrupted trajectory using old/new fixed binaries,
Python mapping and fallback, preparation manifests and sdist→wheel. No fixed
binary edits. New explicit external build, at most4 cores, no remote Actions.
